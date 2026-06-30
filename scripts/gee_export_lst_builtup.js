var grid = ee.FeatureCollection('users/hothikimanh17/hanoi_grid_250m_upload');
var aoi = grid.geometry();
var aoiSimple = aoi.simplify(100);

Map.centerObject(aoi, 10);
Map.addLayer(grid, {}, 'Hanoi grid 250m', false);
print('Grid cell count:', grid.size());

var startDate = '2024-05-01';
var endDate = '2024-08-31';

var inputBands = [
  'B2', 'B3', 'B4',
  'B5', 'B6', 'B7', 'B8', 'B8A',
  'B11', 'B12',
  'NDVI', 'NDWI', 'MNDWI', 'NDBI',
  'NDRE', 'SAVI', 'BSI'
];

function maskS2SCL(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(0)
    .and(scl.neq(1))
    .and(scl.neq(2))
    .and(scl.neq(3))
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));

  return image.updateMask(mask)
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}

// Nâng ngưỡng mây 45% để có đủ ảnh mùa hè (Hà Nội nhiều mây tháng 5-8)
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoiSimple)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 45))
  .map(maskS2SCL);

print('Số ảnh S2 sau bộ lọc mây < 45%:', s2.size());
var composite = s2.median();

var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');
var ndwi = composite.normalizedDifference(['B3', 'B8']).rename('NDWI');
var mndwi = composite.normalizedDifference(['B3', 'B11']).rename('MNDWI');
var ndbi = composite.normalizedDifference(['B11', 'B8']).rename('NDBI');
var ndre = composite.normalizedDifference(['B8A', 'B5']).rename('NDRE');

var savi = composite.expression(
  '((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
  'NIR': composite.select('B8'),
  'RED': composite.select('B4'),
  'L': 0.5
}
).rename('SAVI');

var bsi = composite.expression(
  '((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
  'SWIR1': composite.select('B11'),
  'RED': composite.select('B4'),
  'NIR': composite.select('B8'),
  'BLUE': composite.select('B2')
}
).rename('BSI');

var classifyImage = composite
  .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'])
  .addBands([ndvi, ndwi, mndwi, ndbi, ndre, savi, bsi])
  .rename(inputBands);


var lcValidMask = classifyImage.select('B4').mask().rename('lc_valid');

var coverageCheck = lcValidMask.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: aoiSimple,
  scale: 100,
  maxPixels: 1e9,
  tileScale: 8,
  bestEffort: true
});
print('S2 composite valid coverage May-Aug 2024:', coverageCheck);

// --------------------------------------------------------------------------
// 5. VISUALIZATION
// --------------------------------------------------------------------------

Map.addLayer(
  composite.clip(aoiSimple),
  { bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3 },
  'S2 True Color May-Aug 2024',
  true
);
Map.addLayer(
  composite.clip(aoiSimple),
  { bands: ['B8', 'B4', 'B3'], min: 0, max: 0.4 },
  'S2 False Color NIR-Red-Green',
  false
);
Map.addLayer(ndvi.clip(aoiSimple), { min: -0.2, max: 0.8, palette: ['red', 'yellow', 'green'] }, 'NDVI', false);
Map.addLayer(ndre.clip(aoiSimple), { min: 0, max: 0.6, palette: ['brown', 'yellow', 'green'] }, 'NDRE', false);
Map.addLayer(mndwi.clip(aoiSimple), { min: -0.5, max: 0.5, palette: ['brown', 'white', 'blue'] }, 'MNDWI', false);
Map.addLayer(ndbi.clip(aoiSimple), { min: -0.5, max: 0.5, palette: ['green', 'white', 'red'] }, 'NDBI', false);
Map.addLayer(bsi.clip(aoiSimple), { min: -0.5, max: 0.5, palette: ['green', 'white', 'orange'] }, 'BSI', false);

var waterFC = water.map(function (f) { return f.set('landcover', 0); });
var treeFC = tree.map(function (f) { return f.set('landcover', 1); });
var grassFC = grass.map(function (f) { return f.set('landcover', 2); });
var cropFC = crop.map(function (f) { return f.set('landcover', 3); });
var builtFC = built.map(function (f) { return f.set('landcover', 4); });
var bareFC = bare.map(function (f) { return f.set('landcover', 5); });

var trainingFC = waterFC
  .merge(treeFC)
  .merge(grassFC)
  .merge(cropFC)
  .merge(builtFC)
  .merge(bareFC);

print('Tổng đa giác mẫu đầu vào:', trainingFC.size());

// Chia đa giác ngẫu nhiên 70/30 (seed cố định để tái lập kết quả)
var polygonsWithRandom = trainingFC.randomColumn('poly_random', 42);
var trainPolygons = polygonsWithRandom.filter(ee.Filter.lt('poly_random', 0.7));
var testPolygons = polygonsWithRandom.filter(ee.Filter.gte('poly_random', 0.7));

print('Đa giác train (70%):', trainPolygons.size());
print('Đa giác test  (30%):', testPolygons.size());

// Sample pixel từ mỗi nhóm đa giác riêng biệt
var trainSamplesRaw = classifyImage.sampleRegions({
  collection: trainPolygons,
  properties: ['landcover'],
  scale: 10,
  tileScale: 8,
  geometries: false
});

var testSamplesRaw = classifyImage.sampleRegions({
  collection: testPolygons,
  properties: ['landcover'],
  scale: 10,
  tileScale: 8,
  geometries: false
});

print('Pixel mẫu thô — train:', trainSamplesRaw.size());
print('Pixel mẫu thô — test:', testSamplesRaw.size());
print('Phân bố train thô theo lớp:', trainSamplesRaw.aggregate_histogram('landcover'));
print('Phân bố test  thô theo lớp:', testSamplesRaw.aggregate_histogram('landcover'));

// --------------------------------------------------------------------------
// 7. CÂN BẰNG SỐ MẪU (balanced sampling per class)
// --------------------------------------------------------------------------

function capClassSamples(fc, classId, maxN, seed) {
  return fc
    .filter(ee.Filter.eq('landcover', classId))
    .randomColumn('sample_random_' + classId, seed + classId)
    .sort('sample_random_' + classId)
    .limit(maxN);
}

function buildBalancedSamples(fc, maxN, seed) {
  return capClassSamples(fc, 0, maxN, seed)
    .merge(capClassSamples(fc, 1, maxN, seed))
    .merge(capClassSamples(fc, 2, maxN, seed))
    .merge(capClassSamples(fc, 3, maxN, seed))
    .merge(capClassSamples(fc, 4, maxN, seed))
    .merge(capClassSamples(fc, 5, maxN, seed));
}

// 2000 pixel/lớp train, 700 pixel/lớp test → tổng ~12k + ~4.2k
var trainingSamples = buildBalancedSamples(trainSamplesRaw, 2000, 100);
var testingSamples = buildBalancedSamples(testSamplesRaw, 700, 200);

print('Mẫu train cân bằng (theo lớp):', trainingSamples.aggregate_histogram('landcover'));
print('Mẫu test  cân bằng (theo lớp):', testingSamples.aggregate_histogram('landcover'));
print('Tổng mẫu train:', trainingSamples.size());
print('Tổng mẫu test:', testingSamples.size());

// --------------------------------------------------------------------------
// 8. HUẤN LUYỆN RANDOM FOREST
// --------------------------------------------------------------------------

var classifier = ee.Classifier.smileRandomForest({
  numberOfTrees: 300,
  bagFraction: 0.7,
  seed: 42
}).train({
  features: trainingSamples,
  classProperty: 'landcover',
  inputProperties: inputBands
});

// --------------------------------------------------------------------------
// 9. PHÂN LOẠI & HIỂN THỊ
// --------------------------------------------------------------------------

var classified = classifyImage.classify(classifier).rename('landcover');

var classPalette = [
  '4575b4', // 0 water
  '1a9641', // 1 tree
  'a6d96a', // 2 grass
  'd9ef8b', // 3 crop
  'd73027', // 4 built
  'fdae61'  // 5 bare
];

// Chỉ clip khi hiển thị lên Map, KHÔNG clip image trước khi classify
Map.addLayer(
  classified.clip(aoiSimple),
  { min: 0, max: 5, palette: classPalette },
  'RF Land Cover May-Aug 2024',
  true
);

// --------------------------------------------------------------------------
// 10. ĐÁNH GIÁ ĐỘ CHÍNH XÁC (ACCURACY ASSESSMENT)
// --------------------------------------------------------------------------

var testClassified = testingSamples.classify(classifier);
var confMatrix = testClassified.errorMatrix('landcover', 'classification');

print('=== ACCURACY ASSESSMENT ===');
print('Confusion Matrix:', confMatrix);
print('Overall Accuracy:', confMatrix.accuracy());
print('Kappa coefficient:', confMatrix.kappa());
print('Producer accuracy per class (Recall):', confMatrix.producersAccuracy());
print('User accuracy per class (Precision):', confMatrix.consumersAccuracy());

// --------------------------------------------------------------------------
// 11. VARIABLE IMPORTANCE (RF)
// --------------------------------------------------------------------------

var explain = classifier.explain();
print('RF variable importance:', ee.Dictionary(explain.get('importance')));

// --------------------------------------------------------------------------
// 12. TÍNH PHẦN TRĂM PHÂN LỚP CHO LƯỚI 250M
// --------------------------------------------------------------------------

var waterMask = classified.eq(0).rename('is_water');
var treeMask = classified.eq(1).rename('is_tree');
var grassMask = classified.eq(2).rename('is_grass');
var cropMask = classified.eq(3).rename('is_crop');
var builtMask = classified.eq(4).rename('is_built');
var bareMask = classified.eq(5).rename('is_bare');
var validMask = classifyImage.select('B4').mask().rename('lc_valid');

var lcStack = waterMask
  .addBands(treeMask)
  .addBands(grassMask)
  .addBands(cropMask)
  .addBands(builtMask)
  .addBands(bareMask)
  .addBands(validMask);

var gridLC = lcStack.reduceRegions({
  collection: grid,
  reducer: ee.Reducer.mean(),
  scale: 10,
  tileScale: 8
});

var getVal = function (feat, prop) {
  var v = feat.get(prop);
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), 0, v));
};

var gridFinal = gridLC.map(function (f) {
  var water_pct = getVal(f, 'is_water').multiply(100);
  var tree_pct = getVal(f, 'is_tree').multiply(100);
  var grass_pct = getVal(f, 'is_grass').multiply(100);
  var crop_pct = getVal(f, 'is_crop').multiply(100);
  var built_pct = getVal(f, 'is_built').multiply(100);
  var bare_pct = getVal(f, 'is_bare').multiply(100);
  var lc_valid_pct = getVal(f, 'lc_valid').multiply(100);

  return f.set({
    'water_pct': water_pct,
    'tree_pct': tree_pct,
    'grass_pct': grass_pct,
    'crop_pct': crop_pct,
    'built_pct': built_pct,
    'bare_pct': bare_pct,
    'green_pct': tree_pct.add(grass_pct),
    'physical_greenblue_pct': tree_pct.add(grass_pct).add(water_pct),
    'lc_valid_pct': lc_valid_pct,
    'lc_valid_flag': lc_valid_pct.gte(80)  // 1 = ô lưới có ≥80% pixel hợp lệ
  });
});

print('Grid LC sample (3 ô đầu):', gridFinal.limit(3));

// --------------------------------------------------------------------------
// 13. EXPORTS
// --------------------------------------------------------------------------

// A. Đa giác mẫu huấn luyện (backup / kiểm tra phân bố)
Export.table.toDrive({
  collection: trainingFC,
  description: 'Hanoi_RF_TrainingPolygons_6class_MayAug2024',
  folder: 'GEE_Hanoi_LC',
  fileNamePrefix: 'hanoi_rf_training_polygons_6class_mayaug2024',
  fileFormat: 'SHP'
});

// B. Phần trăm land cover cho từng ô lưới 250m → Step 8 Python
Export.table.toDrive({
  collection: gridFinal.select([
    'grid_id',
    'water_pct', 'tree_pct', 'grass_pct', 'crop_pct',
    'built_pct', 'bare_pct',
    'green_pct', 'physical_greenblue_pct',
    'lc_valid_pct', 'lc_valid_flag',
    '.geo'
  ]),
  description: 'Hanoi_Grid250m_LandCover_MayAug2024_6class_RF_CSV',
  folder: 'GEE_Hanoi_LC',
  // Tên file khớp với config.LANDCOVER_CSV trong Step 8
  fileNamePrefix: 'hanoi_lc_maroct2024_6class_rededge_grid250m',
  fileFormat: 'CSV'
});

// C. Raster phân loại 10m (để kiểm tra thị giác và backup)
Export.image.toDrive({
  image: classified.toByte(),
  description: 'Hanoi_LandCover_MayAug2024_6class_RF_Raster',
  folder: 'GEE_Hanoi_LC',
  fileNamePrefix: 'hanoi_lc_mayaug2024_6class_rf_raster',
  region: aoiSimple,
  scale: 10,
  crs: 'EPSG:32648',
  maxPixels: 1e13
});

// D. Lưới 250m gốc (backup SHP)
Export.table.toDrive({
  collection: grid,
  description: 'Hanoi_Grid250m_Export',
  folder: 'GEE_Hanoi_Grid',
  fileNamePrefix: 'hanoi_grid_250m',
  fileFormat: 'SHP'
});
