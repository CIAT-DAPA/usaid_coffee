

////// Setting the Dates Limits

var lagRange=0
var dates1 =ee.Date('2015-01-01')
var dates2= ee.Date('2016-12-31');

var startMillis = dates1.millis();
var endMillis = dates2.millis();
var lagMillis = lagRange * (1000*60*60*24);

// imagery start date
var startLagMillis = startMillis.subtract(lagMillis);
var startLagDate = ee.Date(startLagMillis);


// imagery end date
var endLagMillis = endMillis.add(lagMillis);
var endLagDate = ee.Date(endLagMillis);


var endMillExclusive = endMillis.add(1000*60*60*24);
print(startLagDate)
print(endLagDate)

/////////////////////////////////////
///////////  FUNCTIONS
////*******
////////////
var maskClouds = function(image) {
  var scored = ee.Algorithms.Landsat.simpleCloudScore(image);
  return image.updateMask(scored.select(['cloud']).lt(30));
};
var cloud_mask = function(image) {
  return image.mask(image.select('BQA').lte(20516))
  .addBands(image.metadata('system:time_start'));
};

var BandsCond = function(image){
  return image.addBands(image.expression(
    '(0.059< swir1Cond && swir1Cond < 0.15) && (0.10< blueCond && blueCond < 0.25)&& (0.45< ndCond && ndCond < 0.6)', {
      'swir1Cond': image.select('swir1'),
      'blueCond': image.select('blue'),
      'ndCond': image.select('nd')
    }).rename('defCond'))
}


var BSI = function(image) {
  return image
    // NDVI
    .addBands(image.expression(
    '((b6 + b4) - (b5+b2))/((b6+b4)+ (b5+b2))', {
      'b6': image.select('B6'),
      'b4': image.select('B4'),
      'b5': image.select('B5'),
      'b2': image.select('B2')
    }).rename('BSI'))
};

var addVI = function(image) {
  return image
    // NDVI
    .addBands(image.normalizedDifference(['nir', 'red']))
    .addBands(image.expression(
    '2.5 * ((NIR - RED) / (NIR + (6 * RED) - (7.5 * BLUE) + 1))', {
      'NIR': image.select('nir'),
      'RED': image.select('red'),
      'BLUE': image.select('blue')
    }).rename('EVI'))
    .addBands(image.expression(
    '((NIR - SWIR1) /(NIR + SWIR1) )', {
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
    }).rename('NDMI'))
    // time in days
    .addBands(image.metadata('system:time_start'));
};


var ConditionalMask = function(image) {
  return image.updateMask(image.select(['defCond']).eq(0));
};

var Fmask = function(image) {
  return image.mask(image.select('fmask').eq(0));
};

var HottestMask = function(image) {
  return image.updateMask(image.select(['Hottest']).lt(0.007));
};

var HottestMeanVis = function(image) {
  return image
    .addBands(image.expression(
    'BLUE- 0.5*RED - 0.08', {
      'RED': image.select('red'),
      'BLUE': image.select('blue')
    }).rename('Hottest'))
    .addBands(image.expression(
    '(GREEN + BLUE + RED)/3', {
      'RED': image.select('red'),
      'BLUE': image.select('blue'),
      'GREEN': image.select('green')
    }).rename('MeanVis'))
    .addBands(image.expression(
    '(BLUE / ND)', {
      'ND': image.select('nd'),
      'BLUE': image.select('blue')
    }).rename('bluend'))
    // time in days
    .addBands(image.metadata('system:time_start'));
};



var whitenessTest = function(image) {
  return image
    .addBands(image.expression(
    '(BLUE- MEANVIS)/MEANVIS', {
      'BLUE': image.select('blue'),
      'MEANVIS': image.select('MeanVis')
    }).abs().rename('bluemv'))
    .addBands(image.expression(
    '(GREEN- MEANVIS)/MEANVIS', {
      'GREEN': image.select('green'),
      'MEANVIS': image.select('MeanVis')
    }).abs().rename('greenmv'))
    .addBands(image.expression(
    '(RED- MEANVIS)/MEANVIS', {
      'RED': image.select('red'),
      'MEANVIS': image.select('MeanVis')
    }).abs().rename('redmv'));
};
var whitenessTestSum = function(image) {
  return image
    .addBands(image.expression(
    '(BLUE+ GREEN + RED)', {
      'BLUE': image.select('bluemv'),
      'GREEN': image.select('greenmv'),
      'RED': image.select('redmv')
    }).abs().rename('WT_Value'));
};

var whitenessTestMask = function(image) {
  return image.updateMask(image.select(['WT_Value']).lt(0.78));
};


var Shadowtest = function(image) {
  return image.addBands(image.select('nir').gt(0.18)
             .focal_min({kernel: kernel, iterations: 2})
             .focal_max({kernel: kernel, iterations: 2})
    .rename('testNIR'));
};

var ShadowTestMask = function(image) {
  return image.updateMask(image.select('testNIR'));
};




/////////////////////////////////////////////
///////
/////// Illumination Correction


var terrain = ee.call('Terrain', ee.Image('USGS/SRTMGL1_003'));

var ic_function = function(image){
  var image2=image
  var solar_zenith = ee.Number(image.get('SUN_ELEVATION')).add(-90).multiply(-1)
  var solar_azimuth = ee.Number(image.get('SUN_AZIMUTH'))
  var solar_zenith_radians = solar_zenith.multiply(Math.PI).divide(180);

  var slope_radians = terrain.select(['slope']).expression("(b('slope')*" + Math.PI + ")/180");
  var aspect = terrain.select(['aspect']);
  
  //slope part of the illumination condition
  var cosZ = solar_zenith_radians.cos();

  var cosS = slope_radians.cos();

  var slope_illumination = cosS.expression("Slope * cosZ ",{
      'Slope': cosS.select('slope'),
      'cosZ': cosZ
    }).rename('b1');

  
  //aspect part of the illumination condition
  var sinZ = solar_zenith_radians.sin();
  var sinS = slope_radians.sin();
  var azimuth_diff_radians = aspect.expression("Aspect - solar_azimuth",{
      'Aspect': aspect.select('aspect'),
      'solar_azimuth': solar_azimuth

    });
  
  azimuth_diff_radians= azimuth_diff_radians.multiply(Math.PI).divide(180)
  
  var cosPhi = azimuth_diff_radians.cos();

  var aspect_illumination = cosPhi.multiply(sinS).multiply(sinZ).select(['aspect'],['b1']);

  var ic = slope_illumination.add(aspect_illumination);
  return image.expression("((image * (cosZ + coeff)) / (ic + coeff)) + offsets", {
    'image': image2.select(DEF_BANDS),
    'ic': ic,
    'cosZ': cosZ,
    'coeff': [5.10309592743233,2.01600885697925,1.936130494609,0.986886169680864,0.548283561041587    ,0.548283561041587],
    'offsets': [0, 0, 0,0,0,0]
  }).set('system:time_start', image2.get('system:time_start'));
}


//////////////*********
///////// ------>  Reading the Feature Collections
//////////////*********

var departament = ee.FeatureCollection('ft:1_wQck1_7lSnvh27_R6ezQxsfk2NlkyHYdkNk8Odc');
var Risaralda =ee.Feature(departament.first()).geometry();

//18xIlNDe5zDDmru6yoqtKT8GPwXFV28YnkOmUwh6a
//1MMmVf5EuzeOvms6rnES9AygVStnQzczqdNQ7k8f6 all points
var ID06 = ee.FeatureCollection('ft:1MMmVf5EuzeOvms6rnES9AygVStnQzczqdNQ7k8f6');
print(ID06)

// selecting a coffee plot
var ID05 = ee.Geometry(ID06.filter(ee.Filter.eq('Lote', 'L_660450184503')).geometry());
print(ID05)

//////////////*********
///////// ------>  Loading the Image Collections
//////////////*********

var LC8_BANDS = ['B2',   'B3',    'B4',  'B5',  'B6',  'B7'];
var LC7_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',  'B7'];
var LC5_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',  'B7'];
var DEF_BANDS = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];


var imageL8 = ee.ImageCollection('LANDSAT/LC8_L1T_TOA_FMASK')
 .filterBounds(ID06)
 .map(Fmask)
 .map(cloud_mask)
 .map(maskClouds)
 .filterDate(startLagDate, endLagDate)
 .map(BSI)

imageL8= imageL8.select(LC8_BANDS,DEF_BANDS)
 .map(ic_function);


var imageL7 = ee.ImageCollection('LANDSAT/LE7_L1T_TOA_FMASK')
.filterBounds(Risaralda)
.map(Fmask)
.map(maskClouds)
.filterDate(startLagDate, endLagDate)
.select(LC7_BANDS,DEF_BANDS)
.map(ic_function);
print(imageL7);

var imageL5 = ee.ImageCollection('LANDSAT/LT5_L1T_TOA_FMASK')
.filterBounds(Risaralda)
.map(Fmask)
.map(maskClouds)
.filterDate(startLagDate, endLagDate)
.select(LC5_BANDS,DEF_BANDS)
.map(ic_function);


var L8_L7Coll=ee.ImageCollection(imageL8.merge(imageL7).merge(imageL5)).map(addVI);


L8_L7Coll = L8_L7Coll.sort("system:time_start")

L8_L7Coll=L8_L7Coll.map(BandsCond).map(ConditionalMask);

var L8_L7Coll3=L8_L7Coll.map(HottestMeanVis).map(HottestMask);


/////////////////// whitenessTest Masking

var L8_L7Coll_WT=L8_L7Coll3.map(whitenessTest).map(whitenessTestSum);

L8_L7Coll_WT=L8_L7Coll_WT.map(whitenessTestMask);

/////////////////// Shadow Masking

var kernel = ee.Kernel.circle({radius:2});

var L8_L7Coll_ST=L8_L7Coll_WT.map(Shadowtest).map(ShadowTestMask);



//////////////*********
///////// ------>  Exporting the points values in CSV
//////////////*********

var imageCFree = L8_L7Coll_ST;
print(imageCFree)

var Tseries = imageCFree.map(function(img) {
  var img2= ee.Image(img);
  var mean = img2.reduceRegions(ID06,'mean',10);
  
  return mean.map(function(f) {
    return f.set('date', ee.Date(img2.get('system:time_start')).format('YYYY-MM-dd'))
  })
}).flatten();

var fields = ["date", "OBJECTID",'blue', 'green', 'red', 'nir', 'swir1', 'swir2'];
var result = ee.FeatureCollection([ee.Feature(null, ee.Dictionary.fromLists(fields, fields))])
    .merge(Tseries)

var resultRed=result.filter(ee.Filter.neq('blue', null))

//print (resultRed)
Export.table.toDrive({
  collection: resultRed,
  description:'mean_bands_coffe_IC_93_5mil',
  fileFormat: 'CSV'
});



//////////////*********
///////// ------>  Time Series Visualization
//////////////*********

//Map.centerObject(ID05,16);
var chart = ui.Chart.image.series({
  imageCollection:  L8_L7Coll_ST.select(['system:time_start','nd','swir1','blue','swir2']),
  region: ID05,
  reducer: ee.Reducer.mean(),
  scale: 5
});

chart.style().set({
  position: 'bottom-right',
  width: '800px',
  height: '300px'
});
Map.add(chart);

var sfLayer = ui.Map.Layer(ID05, {color: 'FF0000'}, 'Coffe');
Map.layers().add(sfLayer);



// Create a label on the map.
var label = ui.Label('Click a point on the chart to show the image for that date.');
Map.add(label);

// When the chart is clicked, update the map and label.

chart.onClick(function(xValue, yValue, seriesName) {
  if (!xValue) return;  // Selection was cleared.

  // Show the image for the clicked date.
  var equalDate = ee.Filter.equals('system:time_start', xValue);
  var image = ee.Image(L8_L7Coll_ST.filter(equalDate).first());
  var l8Layer = ui.Map.Layer(image, {
    gamma: 1.3,
    min: 0,
    max: 0.3,
    bands: ['red', 'green', 'blue']
  });
  var image2 = ee.Image(L8_L7Coll.filter(equalDate).first());
  var l8Layer2 = ui.Map.Layer(image2, {
    gamma: 1.3,
    min: 0,
    max: 0.3,
    bands: ['red', 'green', 'blue']
  });
  Map.layers().reset([l8Layer2,l8Layer, sfLayer]);
  // Show a label with the date on the map.
  label.setValue((new Date(xValue)).toUTCString());
  print(xValue);
});


///////////// Validated
var medianLandsat7 = imageL7.median();
var medianLandsat8 = imageL8.median();
//Add to map

var vizParams = {bands: ['red','green','blue'], min: 0, max: 0.4};
Map.addLayer(medianLandsat7, vizParams, 'medianLandsat7', 0);
Map.addLayer(medianLandsat8, vizParams, 'medianLandsat8', 0);

//Extract reflectance values over random points for L7 and L8 bands
var join = medianLandsat7.addBands(medianLandsat8);
print(join, 'joint image');

var proj = ee.Image(imageL7.first()).projection();

Map.addLayer(ID05.buffer(4000), {}, 'area');
var pointNum = ee.Number(100);//Adjust this number depending on the desired number of points
var rdm = ee.FeatureCollection.randomPoints(ID05.buffer(4000), pointNum);
Map.addLayer(rdm, {}, 'random points', 0);


var extract = join.reduceRegions({
  collection: rdm,
  reducer: ee.Reducer.mean(),
  scale: 20,
});

print(extract.limit(20), 'random points extraction');

// Export if needed
//Export.table(ee.FeatureCollection(extract).flatten());

// Create chart to visualize relationship
var chart = ui.Chart.feature.byFeature(extract,'blue', 'blue_1').setChartType('ScatterChart')
    .setOptions({title: 'S2 vs S1',
    hAxis: {'title': 'landsat 7 reflectance'},
    vAxis: {'title': 'landsat 8 reflectance'},
    pointSize: 1,});
print(chart);


