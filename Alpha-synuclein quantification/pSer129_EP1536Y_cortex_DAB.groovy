
// Script developed by Irene Frigerio (2021)
// Used in ex. Frigerio et al. 2023, Translational Neurodegeneration PMID: 36658627
// QuPath version 0.2.3
// Microscope: whole-slide scanner Vectra Polaris, 20x  objective
// Antibody pSer129 aSyn, clone EP153Y6 (Abcam, article nr ab51253)
// Please do not forget to adjust the output folder


// This script quantifies Lewy Bodies (LBs) and remaining alpha-synuclein pathology (% area) in the cortex
// Strategy for this script: 
//    1.  Detect Lewy bodies
//    2.  Area% of alpha-synuclein after excluding LBs 

// Standard settings, set deconvolutions for DAB and hematoxylin
setImageType('BRIGHTFIELD_H_DAB');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759 ", "Background" : " 255 255 255 "}');


                // 1 QUANTIFY LEWY BODY COUNT (as outcome measure normalize it over mm2)
// 1.1 Detect Lewy- like body-like objects based on optimal density sum
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield": "Optical density sum",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 25.0,  "medianRadiusMicrons": 1.0,  "sigmaMicrons": 3.0,  "minAreaMicrons": 19.6,  "maxAreaMicrons": 490.9,  "threshold": 0.6,  "maxBackground": 8.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

// 1.2 Select only objects that are DAB positive, delete the other objects
runObjectClassifier("LB_DAB");
selectObjectsByClassification("Background");
clearSelectedObjects(true);
selectObjects();

// 1.3 Count lewy bodies  ---- **** TO SAVE YOUR OWN DATA, CORRECT THIS FILE LOCATION
saveAnnotationMeasurements('/C:/path/to/output/folder/LBs')


// 1.4 Make from the object detections (lewy bodies) annotations 
def detectionsLewy = getDetectionObjects().findAll{it.getPathClass() == getPathClass("Lewy Body")}
def newAnnotationsLewy = detectionsLewy.collect {
    return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass()) 
}
removeObjects(detectionsLewy, true)
insertObjects(newAnnotationsLewy)

print "done with Lewy Bodies" 

            // 2. QUANTIFY REMAINING P-SER129 PATHOLOGY (% AREA)
// 2.1 Substract LB annotations from ROI
resolveHierarchy();

import qupath.lib.roi.*
import qupath.lib.objects.*

classToSubtract = 'Lewy Body'
    
def topLevel = getObjects{return it.getLevel()==1 && it.isAnnotation()}
println(topLevel)
for (parent in topLevel){

    def total = []
    def polygons = []
    subtractions = parent.getChildObjects().findAll{it.isAnnotation() }
    println(subtractions)
    for (subtractyBit in subtractions){
        if (subtractyBit instanceof AreaROI){
           subtractionROIs =splitAreaToPolygons(subtractyBit.getROI())
           total.addAll(subtractionROIs[1])
        } else {total.addAll(subtractyBit.getROI())}              
                
    }     
    if (parent instanceof AreaROI){
        polygons = RoiTools.splitAreaToPolygons(parent.getROI())
        total.addAll(polygons[0])
    } else { polygons[1] = parent.getROI()}

            
    def newPolygons = polygons[1].collect {
    updated = it
    for (hole in total)
         updated = RoiTools.combineROIs(updated, hole, RoiTools.CombineOp.SUBTRACT)
         return updated
    }
                // Remove original annotation, add new ones
    annotations = newPolygons.collect {new PathAnnotationObject(updated, parent.getPathClass())}


    addObjects(annotations)

    removeObjects(subtractions, true)
    removeObject(parent, true)
}
print "Done with substracting Lewy Bodies from ROI"

// 2.2 pSer129 aSyn pixel classifier 
selectAnnotations();
addPixelClassifierMeasurements("aSyn_pathology", "aSyn_pathology")
saveAnnotationMeasurements('/C:/path/to/output/folder/non-LB')

print "Done with script"