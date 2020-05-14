// Javascript file to extract environmental data from Crowther Lab Composite on Google Earth Engine
// Covariates to sample
var composite = ee.Image('users/crowtherlab/Composite/CrowtherLab_Composite_30ArcSec'); 
var resolveBiomes = ee.Image("users/devinrouth/Resolve_Biomes_30ArcSec");

// Load featurecollection
var data_fc = ee.FeatureCollection("users/johanvandenhoogen/colin_temp/data_for_composite_CA")
Map.addLayer(data_fc, {}, 'Sampling points')
print('Size of FC', data_fc.size())

var compositeToSample = composite.addBands(resolveBiomes.toInt());

var sampledPointsToExport = data_fc.map(function(f){
  return f.set(compositeToSample.reduceRegion('first',f.geometry()));
});

Export.table.toDrive({
	collection: sampledPointsToExport,
	description: '20200505_CA_dataset_Sampled'
});