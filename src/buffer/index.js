import clone from '../clone';
import { BufferOp, GeoJSONReader, GeoJSONWriter } from 'turf-jsts';
import centerOfMass from '../center-of-mass';
import { geomEach, coordEach, featureEach } from '../meta';
import { featureCollection,
    radiansToLength,
    lengthToRadians,
    feature,
    degreesToRadians, radiansToDegrees, checkIfOptionsExist } from '../helpers';

/**
 * Calculates a buffer for input features for a given radius. Units supported are miles, kilometers, and degrees.
 *
 * When using a negative radius, the resulting geometry may be invalid if
 * it's too small compared to the radius magnitude. If the input is a
 * FeatureCollection, only valid members will be returned in the output
 * FeatureCollection - i.e., the output collection may have fewer members than
 * the input, or even be empty.
 *
 * @name buffer
 * @param {FeatureCollection|Geometry|Feature<any>} geojson input to be buffered
 * @param {number} radius distance to draw the buffer (negative values are allowed)
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units="kilometers"] any of the options supported by turf units
 * @returns {FeatureCollection|Feature<Polygon|MultiPolygon>|undefined} buffered features
 * @example
 * var point = turf.point([-90.548630, 14.616599]);
 * var buffered = turf.buffer(point, 500, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [point, buffered]
 */
function buffer(geojson, radius, options) {

    // Optional params
    options = checkIfOptionsExist(options);
    var units = options.units || 'kilometers';
    // var steps = options.steps || 64;

    // validation
    if (!geojson) throw new Error('geojson is required');
    if (typeof options !== 'object') throw new Error('options must be an object');
    // if (typeof steps !== 'number') throw new Error('steps must be an number');

    // Allow negative buffers ("erosion") or zero-sized buffers ("repair geometry")
    if (radius === undefined) throw new Error('radius is required');
    // if (steps <= 0) throw new Error('steps must be greater than 0');

    var distance = radiansToLength(lengthToRadians(radius, units), 'meters');
    var results = [];

    switch (geojson.type) {
        case 'GeometryCollection':
            geomEach(geojson, function (geometry) {
                var buffered = bufferFeature(geometry, distance);
                if (buffered) results.push(buffered);
            });
            return featureCollection(results);
        case 'FeatureCollection':
            featureEach(geojson, function (feature) {
                var buffered = bufferFeature(feature, distance);
                if (buffered) results.push(buffered);
            });
            return featureCollection(results);
    }
    return bufferFeature(geojson, distance);
}

/**
 * Buffer single Feature/Geometry
 *
 * @private
 * @param {Feature<any>} geojson input to be buffered
 * @param {number} radius distance to draw the buffer
 * @param {number} [steps=64] number of steps
 * @returns {Feature<Polygon|MultiPolygon>} buffered feature
 */
function bufferFeature(geojson, radius) {
    var properties = geojson.properties || {};
    var geometry = (geojson.type === 'Feature') ? clone(geojson.geometry) : clone(geojson);

    var centroid = centerOfMass(geometry);
    var utmZone = checkUtmZone(centroid.geometry);
    reprojectFeature(geometry, utmZone, true);

    var reader = new GeoJSONReader();
    var geom = reader.read(geometry);
    var buffered = BufferOp.bufferOp(geom, radius);
    var writer = new GeoJSONWriter();
    buffered = writer.write(buffered);

    if (coordsIsNaN(buffered.coordinates)) return undefined;
    reprojectFeature(buffered, utmZone, false);

    return (buffered.geometry) ? buffered : feature(buffered, properties);
}

export default buffer;

/**
 * Coordinates isNaN
 *
 * @private
 * @param {Array<any>} coords GeoJSON Coordinates
 * @returns {boolean} if NaN exists
 */
function coordsIsNaN(coords) {
    if (Array.isArray(coords[0])) return coordsIsNaN(coords[0]);
    return isNaN(coords[0]);
}

function reprojectFeature(feature, utmZone, toUtm) {
    coordEach(feature, function (coord, coordIndex) { //eslint-disable-line
        var blah = toUtm ? convertCoordToUtm(coord[0], coord[1], utmZone) : convertUtmToLatLon(coord[0], coord[1], utmZone);
        coord.length = 0;
        coord.push(blah[0], blah[1]);
    }, false);
}

function checkUtmZone(centerPoint) {

    const lat = centerPoint.coordinates[1];
    const lon = centerPoint.coordinates[0];
    let zoneNumber = Math.floor((lon + 180) / 6) + 1;
    let hemisphere = 'N';

    if (lon === 180) zoneNumber = 60;
    if (lat < 0.0) hemisphere = 'S';

    if (lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0) zoneNumber = 32;

    return {zoneNumber, hemisphere};
}

function convertCoordToUtm(lon, lat, zone) {

    let falseEasting = 500e3;
    let falseNorthing = 10000e3;
    let lambda0 = degreesToRadians(((zone.zoneNumber - 1) * 6 - 180 + 3));

    let mgrsLatBands = 'CDEFGHJKLMNPQRSTUVWXX'; // X is repeated for 80-84°N
    let latBand = mgrsLatBands.charAt(Math.floor(lat / 8 + 10));

    // adjust zone & central meridian for Norway
    if (zone === 31 && latBand === 'V' && lon >= 3) { zone++; degreesToRadians(lambda0 += 6); }
    // adjust zone & central meridian for Svalbard
    if (zone === 32 && latBand === 'X' && lon <  9) { zone--; degreesToRadians(lambda0 -= 6); }
    if (zone === 32 && latBand === 'X' && lon >= 9) { zone++; degreesToRadians(lambda0 += 6); }
    if (zone === 34 && latBand === 'X' && lon < 21) { zone--; degreesToRadians(lambda0 -= 6); }
    if (zone === 34 && latBand === 'X' && lon >= 21) { zone++; degreesToRadians(lambda0 += 6); }
    if (zone === 36 && latBand === 'X' && lon < 33) { zone--; degreesToRadians(lambda0 -= 6); }
    if (zone === 36 && latBand === 'X' && lon >= 33) { zone++; degreesToRadians(lambda0 += 6); }

    var fai = degreesToRadians(lat);      // latitude ± from equator
    var lambda = degreesToRadians(lon) - lambda0; // longitude ± from central meridian

    let a = 6378137;
    let f = 1 / 298.257223563;
    // WGS 84: a = 6378137, b = 6356752.314245, f = 1/298.257223563;

    let k0 = 0.9996; // UTM scale on the central meridian

    // ---- easting, northing: Karney 2011 Eq 7-14, 29, 35:

    let e = Math.sqrt(f * (2 - f)); // eccentricity
    let n = f / (2 - f);        // 3rd flattening
    let n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5; // TODO: compare Horner-form accuracy?

    let coslambda = Math.cos(lambda), sinlambda = Math.sin(lambda);

    let tao = Math.tan(fai); // tao ≡ tanfai, tao_ ≡ tanfai_; prime (_) indicates angles on the conformal sphere
    let sigma = Math.sinh(e * Math.atanh(e * tao / Math.sqrt(1 + tao * tao)));

    let tao_ = tao * Math.sqrt(1 + sigma * sigma) - sigma * Math.sqrt(1 + tao * tao);

    let ksai_ = Math.atan2(tao_, coslambda);
    let yita_ = Math.asinh(sinlambda / Math.sqrt(tao_ * tao_ + coslambda * coslambda));

    let A = a / (1 + n) * (1 + 1 / 4 * n2 + 1 / 64 * n4 + 1 / 256 * n6); // 2πA is the circumference of a meridian

    let alpha = [null, // note alpha is one-based array (6th order Krüger expressions)
        1 / 2 * n - 2 / 3 * n2 + 5 / 16 * n3 +   41 / 180 * n4 -     127 / 288 * n5 +      7891 / 37800 * n6,
        13 / 48 * n2 -  3 / 5 * n3 + 557 / 1440 * n4 +     281 / 630 * n5 - 1983433 / 1935360 * n6,
        61 / 240 * n3 -  103 / 140 * n4 + 15061 / 26880 * n5 +   167603 / 181440 * n6,
        49561 / 161280 * n4 -     179 / 168 * n5 + 6601661 / 7257600 * n6,
        34729 / 80640 * n5 - 3418889 / 1995840 * n6,
        212378941 / 319334400 * n6];

    let ksai = ksai_;
    for (let j = 1; j <= 6; j++) ksai += alpha[j] * Math.sin(2 * j * ksai_) * Math.cosh(2 * j * yita_);

    let yita = yita_;
    for (let j = 1; j <= 6; j++) yita += alpha[j] * Math.cos(2 * j * ksai_) * Math.sinh(2 * j * yita_);

    let x = k0 * A * yita;
    let y = k0 * A * ksai;

    x = x + falseEasting;
    if (y < 0) y = y + falseNorthing;

    return [y, x];
}

function convertUtmToLatLon(y, x, zone) {
    var z = zone.zoneNumber;
    var h = zone.hemisphere;

    var falseEasting = 500e3, falseNorthing = 10000e3;

    var a = 6378137, f = 1 / 298.257223563;

    var k0 = 0.9996; // UTM scale on the central meridian

    x = x - falseEasting;               // make x ± relative to central meridian
    y = h === 'S' ? y - falseNorthing : y; // make y ± relative to equator

    // ---- from Karney 2011 Eq 15-22, 36:

    var e = Math.sqrt(f * (2 - f)); // eccentricity
    var n = f / (2 - f);        // 3rd flattening
    var n2 = n * n, n3 = n * n2, n4 = n * n3, n5 = n * n4, n6 = n * n5;

    var A = a / (1 + n) * (1 + 1 / 4 * n2 + 1 / 64 * n4 + 1 / 256 * n6); // 2πA is the circumference of a meridian

    var yita = x / (k0 * A);
    var ksai = y / (k0 * A);

    var beta = [null, // note beta is one-based array (6th order Krüger expressions)
        1 / 2 * n - 2 / 3 * n2 + 37 / 96 * n3 -    1 / 360 * n4 -   81 / 512 * n5 +    96199 / 604800 * n6,
        1 / 48 * n2 +  1 / 15 * n3 - 437 / 1440 * n4 +   46 / 105 * n5 - 1118711 / 3870720 * n6,
        17 / 480 * n3 -   37 / 840 * n4 - 209 / 4480 * n5 +      5569 / 90720 * n6,
        4397 / 161280 * n4 -   11 / 504 * n5 -  830251 / 7257600 * n6,
        4583 / 161280 * n5 -  108847 / 3991680 * n6,
        20648693 / 638668800 * n6];

    var ksai_ = ksai;
    for (var j = 1; j <= 6; j++) ksai_ -= beta[j] * Math.sin(2 * j * ksai) * Math.cosh(2 * j * yita);

    var yita_ = yita;
    for (var j = 1; j <= 6; j++) yita_ -= beta[j] * Math.cos(2 * j * ksai) * Math.sinh(2 * j * yita);  //eslint-disable-line

    var sinhyita_ = Math.sinh(yita_);
    var sinksai_ = Math.sin(ksai_), cosksai_ = Math.cos(ksai_);

    var tao_ = sinksai_ / Math.sqrt(sinhyita_ * sinhyita_ + cosksai_ * cosksai_);

    var taoi = tao_;
    do {
        var sigmai = Math.sinh(e * Math.atanh(e * taoi / Math.sqrt(1 + taoi * taoi)));
        var taoi_ = taoi * Math.sqrt(1 + sigmai * sigmai) - sigmai * Math.sqrt(1 + taoi * taoi);
        var deltataoi = (tao_ - taoi_) / Math.sqrt(1 + taoi_ * taoi_) *
            (1 + (1 - e * e) * taoi * taoi) / ((1 - e * e) * Math.sqrt(1 + taoi * taoi));
        taoi += deltataoi;
    } while (Math.abs(deltataoi) > 1e-12); // using IEEE 754 deltataoi -> 0 after 2-3 iterations
    // note relatively large convergence test as deltataoi toggles on ±1.12e-16 for eg 31 N 400000 5000000
    var tao = taoi;

    var fai = Math.atan(tao);

    var lambda = Math.atan2(sinhyita_, cosksai_);

    var lambda0 = degreesToRadians((z - 1) * 6 - 180 + 3); // longitude of central meridian
    lambda += lambda0; // move lambda from zonal to global coordinates

    // round to reasonable precision
    var lat = radiansToDegrees(fai); // nm precision (1nm = 10^-11°)
    var lon = radiansToDegrees(lambda); // (strictly lat rounding should be fai⋅cosfai!)

    return [lon, lat];
}