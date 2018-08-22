import { flattenEach, coordEach } from '../meta';
import { getCoords, getType } from '../invariant';
import { isObject, lineString, multiLineString, convertLength, degreesToRadians, radiansToDegrees } from '../helpers';
import intersection from './lib/intersection';
import centerOfMass from '../center-of-mass';

/**
 * Takes a {@link LineString|line} and returns a {@link LineString|line} at offset by the specified distance.
 *
 * @name lineOffset
 * @param {Geometry|Feature<LineString|MultiLineString>} geojson input GeoJSON
 * @param {number} distance distance to offset the line (can be of negative value)
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, kilometers, inches, yards, meters
 * @returns {Feature<LineString|MultiLineString>} Line offset from the input line
 * @example
 * var line = turf.lineString([[-83, 30], [-84, 36], [-78, 41]], { "stroke": "#F00" });
 *
 * var offsetLine = turf.lineOffset(line, 2, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [offsetLine, line]
 * offsetLine.properties.stroke = "#00F"
 */
function lineOffset(geojson, distance, options) {
    // Optional parameters
    options = options || {};
    if (!isObject(options)) throw new Error('options is invalid');
    var units = options.units ? options.units : 'kilometers';

    // Valdiation
    if (!geojson) throw new Error('geojson is required');
    if (distance === undefined || distance === null || isNaN(distance)) throw new Error('distance is required');
    const distanceMeters = convertLength(distance, units, 'meters');

    var type = getType(geojson);
    var properties = geojson.properties;

    geojson = JSON.parse(JSON.stringify(geojson));

    switch (type) {
    case 'LineString':
        return lineOffsetFeature(geojson, distanceMeters);
    case 'MultiLineString':
        var coords = [];
        flattenEach(geojson, function (feature) {
            coords.push(lineOffsetFeature(feature, distanceMeters).geometry.coordinates);
        });
        return multiLineString(coords, properties);
    default:
        throw new Error('geometry ' + type + ' is not supported');
    }
}

/**
 * Line Offset
 *
 * @private
 * @param {Geometry|Feature<LineString>} line input line
 * @param {number} distance distance to offset the line (can be of negative value)
 * @param {string} [units=kilometers] units
 * @returns {Feature<LineString>} Line offset from the input line
 */
function lineOffsetFeature(line, distance) {

    var centroid = centerOfMass(line);
    var utmZone = checkUtmZone(centroid.geometry);
    reprojectFeature(line, utmZone, true);

    var segments = [];
    var coords = getCoords(line);
    var finalCoords = [];
    coords.forEach(function (currentCoords, index) {
        if (index !== coords.length - 1) {
            var segment = processSegment(currentCoords, coords[index + 1], distance);
            segments.push(segment);
            if (index > 0) {
                var seg2Coords = segments[index - 1];
                var intersects = intersection(segment, seg2Coords);

                // Handling for line segments that aren't straight
                if (intersects !== false) {
                    seg2Coords[1] = intersects;
                    segment[0] = intersects;
                }

                finalCoords.push(seg2Coords[0]);
                if (index === coords.length - 2) {
                    finalCoords.push(segment[0]);
                    finalCoords.push(segment[1]);
                }
            }
            // Handling for lines that only have 1 segment
            if (coords.length === 2) {
                finalCoords.push(segment[0]);
                finalCoords.push(segment[1]);
            }
        }
    });
    var out = lineString(finalCoords, line.properties);
    reprojectFeature(out, utmZone, false);
    return out;
}

/**
 * Process Segment
 * Inspiration taken from http://stackoverflow.com/questions/2825412/draw-a-parallel-line
 *
 * @private
 * @param {Array<number>} point1 Point coordinates
 * @param {Array<number>} point2 Point coordinates
 * @param {number} offset Offset
 * @returns {Array<Array<number>>} offset points
 */
function processSegment(point1, point2, offset) {
    var L = Math.sqrt((point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]));

    var out1x = point1[0] + offset * (point2[1] - point1[1]) / L;
    var out2x = point2[0] + offset * (point2[1] - point1[1]) / L;
    var out1y = point1[1] + offset * (point1[0] - point2[0]) / L;
    var out2y = point2[1] + offset * (point1[0] - point2[0]) / L;
    return [[out1x, out1y], [out2x, out2y]];
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
                                                             20648693 / 638668800 * n6 ];

    var ksai_ = ksai;
    for (var j = 1; j <= 6; j++) ksai_ -= beta[j] * Math.sin(2 * j * ksai) * Math.cosh(2 * j * yita);

    var yita_ = yita;
    for (var j = 1; j <= 6; j++) yita_ -= beta[j] * Math.cos(2 * j * ksai) * Math.sinh(2 * j * yita);

    var sinhyita_ = Math.sinh(yita_);
    var sinksai_ = Math.sin(ksai_), cosksai_ = Math.cos(ksai_);

    var tao_ = sinksai_ / Math.sqrt(sinhyita_ * sinhyita_ + cosksai_ * cosksai_);

    var taoi = tao_;
    do {
        var sigmai = Math.sinh(e * Math.atanh(e * taoi / Math.sqrt(1 + taoi * taoi)));
        var taoi_ = taoi * Math.sqrt(1 + sigmai * sigmai) - sigmai * Math.sqrt(1 + taoi * taoi);
        var deltataoi = (tao_ - taoi_) / Math.sqrt(1 + taoi_ * taoi_)
            * (1 + (1 - e * e) * taoi * taoi) / ((1 - e * e) * Math.sqrt(1 + taoi * taoi));
        taoi += deltataoi;
    } while (Math.abs(deltataoi) > 1e-12); // using IEEE 754 deltataoi -> 0 after 2-3 iterations
    // note relatively large convergence test as deltataoi toggles on ±1.12e-16 for eg 31 N 400000 5000000
    var tao = taoi;

    var fai = Math.atan(tao);

    var lambda = Math.atan2(sinhyita_, cosksai_);

    // ---- convergence: Karney 2011 Eq 26, 27

    var p = 1;
    for (var j = 1; j <= 6; j++) p -= 2 * j * beta[j] * Math.cos(2 * j * ksai) * Math.cosh(2 * j * yita);
    var q = 0;
    for (var j = 1; j <= 6; j++) q += 2 * j * beta[j] * Math.sin(2 * j * ksai) * Math.sinh(2 * j * yita);

    var gamma_ = Math.atan(Math.tan(ksai_) * Math.tanh(yita_));
    var gamma__ = Math.atan2(q, p);

    var gamma = gamma_ + gamma__;

    // ---- scale: Karney 2011 Eq 28

    var sinfai = Math.sin(fai);
    var k_ = Math.sqrt(1 - e * e * sinfai * sinfai) * Math.sqrt(1 + tao * tao) * Math.sqrt(sinhyita_ * sinhyita_ + cosksai_ * cosksai_);
    var k__ = A / a / Math.sqrt(p * p + q * q);

    var k = k0 * k_ * k__;

    // ------------

    var lambda0 = degreesToRadians((z - 1) * 6 - 180 + 3); // longitude of central meridian
    lambda += lambda0; // move lambda from zonal to global coordinates

    // round to reasonable precision
    var lat = radiansToDegrees(fai); // nm precision (1nm = 10^-11°)
    var lon = radiansToDegrees(lambda); // (strictly lat rounding should be fai⋅cosfai!)

    return [lon, lat];
}

export default lineOffset;
