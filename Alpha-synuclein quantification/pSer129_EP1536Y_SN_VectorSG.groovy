
// Script developed by Chen-Pei Lin (2021)
// Used in Lin et al., 2023, Movement Disorders PMID: 37347552
// QuPath version 0.2.3
// Microscope: whole-slide scanner Olympus VS200 Evident, 20x objective
// Antibody: pSer129 aSyn, clone EP153Y6 (Abcam, article nr ab51253), visualized with Vector SG grey (Vector, California) followed by counter-staining with nuclear fast red (Vector, California)
// Please do not forget to adjust the output folder



// This script quantifies Lewy Bodies (LBs) and rest of alpha-synuclein pathology (% area) in the substantia nigra
// Strategy for this script: 
//    1.  Detect neurons
//    2.  Detect Lewy bodies
//    3.  Include only intracellular lewy bodies
//    4.  Run a Pixel classifier over the ROIs 


// ********************************NECESARRY BEFORE RUNNING THE SCRIPT!! : Classify all ROIs you have drawn to 'parent' **********************************************
// otherwise your annotations will be deleted 


setImageType('BRIGHTFIELD_H_E');
setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "Eosin", "Values 2" : "0.2159 0.8012 0.5581 ", "Background" : " 255 255 255 "}');


   // ***************2. Detect Lewy bodies***************
selectObjectsByClassification("Parent");
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 1.0,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 4.0,  "sigmaMicrons": 3.5,  "minAreaMicrons": 28.27,  "maxAreaMicrons": 99999.0,  "threshold": 0.2,  "maxBackground": 5.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 1.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
runObjectClassifier("Nucleus min caliper");
selectObjectsByClassification("Extracellular Lewy Body");
selectObjectsByClassification("Background");
clearSelectedObjects(true);
clearSelectedObjects();
runObjectClassifier("Nucleus Hematoxylin OD min");
selectObjectsByClassification("Background");
clearSelectedObjects(true);
clearSelectedObjects();

saveAnnotationMeasurements('/C:/path/to/output/folder/LB count/')


fireHierarchyUpdate()
// ************** 4.  Run a Pixel classifier over the ROIs  (Analyse asyn Pathology load)**************

def detectionsLewy = getDetectionObjects().findAll{it.getPathClass() == getPathClass("Lewy Body")}
def newAnnotationsLewy = detectionsLewy.collect {
    return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass()) 
}
removeObjects(detectionsLewy, true)
insertObjects(newAnnotationsLewy)



print "done with Intracellular Lewy Bodies" 

          
// Substract annotations from ROI
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
print "Done with substracting intracellular Lewy Bodies from ROI"

selectAnnotations();
addPixelClassifierMeasurements("aSyn_pathology", "aSyn_pathology")
saveAnnotationMeasurements('/C:/path/to/output/folder/aSyn pathology/')
 
 print "Done with script" 

