//This script quantifies TH neurons (count) and threads (% area), and neuromelanin neurons (count) in tyrosine hydroxilase (TH) staining with vector SG in the substantia nigra

// Script developed by Chen-Pei Lin and Lydian Knoop (2021)
// Used in Lin et al. 2023, Movement Disorders PMID: 37347552
// QuPath version 0.2.3
// Microscope: whole-slide scanner Olympus VS200 Evident, 20x objective
// Antibody TH (Immunostar, article nr 22941)
// Please do not forget to adjust the output folder


//Standard settings. TH neurons in the SN are stained with Vector SG gray and Fast Red, to avoid a similar staining color as neuromelanin 
setImageType('BRIGHTFIELD_OTHER');
setColorDeconvolutionStains('{"Name" : "Vector SG grey - Fast red", "Stain 1" : "Vector SG grey", "Values 1" : "0.61528 0.61528 0.49282 ", "Stain 2" : "Fast Red", "Values 2" : "0.24276 0.71226 0.6586 ", "Background" : " 255 255 255 "}');


// Find Neurons in SN. - With or Without neuromelanin, But with TH for sure
print "Finding TH positive neurons ..." 
selectObjectsByClassification("Parent");
// Find all nuclei
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 1.0,  "backgroundRadiusMicrons": 20.0,  "medianRadiusMicrons": 6.0,  "sigmaMicrons": 3.0,  "minAreaMicrons": 100.0,  "maxAreaMicrons": 9999.0,  "threshold": 0.05,  "maxBackground": 5.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 3.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
runObjectClassifier("Area");
selectObjectsByClassification("Background");
clearSelectedObjects(true);
clearSelectedObjects();
// deletes the threads mistakently picked up before
runObjectClassifier("Nucleus Circularity");
selectObjectsByClassification(null);
clearSelectedObjects(true);
clearSelectedObjects();
runObjectClassifier("Nucleus Fast Red OD range");
selectObjectsByClassification(null);
clearSelectedObjects(true);
clearSelectedObjects();
runObjectClassifier("Cell Fast Red OD std dev");
selectObjectsByClassification(null);
clearSelectedObjects(true);
clearSelectedObjects();
print "...Found them!"

//Save the TH positive neurons 
saveAnnotationMeasurements('/C:/path/to/output/folder/TH_neurons/')

//Make from 'TH neurons' annotations

print "Processing TH neurons..." 
def detectionsTH = getDetectionObjects().findAll{it.getPathClass() == getPathClass("Neurons")}
def newAnnotationsTH = detectionsTH.collect {
    return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass()) 
}
removeObjects(detectionsTH, true)
insertObjects(newAnnotationsTH)

def annotations = getAnnotationObjects()
removeObjects(annotations, false)
addObjects(annotations)


//Find Neuromelanin cells

print "Find Neuromelanin..."
selectObjectsByClassification("Parent");
createDetectionsFromPixelClassifier("Neuromelanin detection", 90.0, 0.0, "SPLIT")
selectObjectsByClassification("Background");
clearSelectedObjects(true);
clearSelectedObjects();

selectObjectsByClassification("Neuromelanin");
def detectionsNM = getDetectionObjects().findAll{it.getPathClass() == getPathClass("Neuromelanin")}


//Remove Neuromalanin cells with TH, so only Neuromenalin remain
print "...Remove the neuromelanins with TH, we don't want them..."
import static qupath.lib.gui.scripting.QPEx.*
def islets=getDetectionObjects()
selectObjectsByClassification("Neurons")
mergeSelectedAnnotations()

def outsidegeo=getAnnotationObjects().find{it.getPathClass()==getPathClass("Neurons")}.getROI().getGeometry()

def intersections=[]
islets.eachWithIndex{entry,idx->
entrygeo = entry.getROI().getGeometry()
    if (!entrygeo.intersects(outsidegeo)){
        intersections<<idx
    }
}

print(intersections)

getCurrentHierarchy().getSelectionModel().selectObjects(islets[intersections], true)
removeObjects(islets[intersections], true)

fireHierarchyUpdate();

print "... Done with Deleting Neuromelanin TH positive neurons"

//Save Neuromelanin only positive neurons 
saveAnnotationMeasurements('/C:/path/to/output/folder/Neuromelanin only neurons/')

//Make from 'Neuromelanin neurons' annotations
print "Starting with extracting all neurons from ROI"
def detectionsNM = getDetectionObjects().findAll{it.getPathClass() == getPathClass("Neuromelanin")}
def newAnnotationsNM = detectionsNM.collect {
    return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass()) 
}
removeObjects(detectionsNM, true)
insertObjects(newAnnotationsNM)

def annotationsNM1 = getAnnotationObjects()
removeObjects(annotationsNM1, false)
addObjects(annotationsNM1)

//Merge both annotations as one, because we want to exclude them in our next step
selectObjectsByClassification("Neurons", "Neuromelanin");
mergeSelectedAnnotations()

// Name the annotations neurons, otherwise we cannot substract them from the ROI
var annotationsz = getAnnotationObjects().findAll{it.getPathClass() == null}
print ("Number of unclassified annotations: " + annotationsz.size())
def classification = getPathClass ('Neuromelanin')
annotationsz.each{it.setPathClass(classification)}

fireHierarchyUpdate()


         
// Substract annotations from ROI
selectObjectsByClassification("Neuromelanin", "Neuromelanin");
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}');

resolveHierarchy();

import qupath.lib.roi.* 
import qupath.lib.objects.*

classToSubtract = 'Neuromelanin'
    
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
print "Done with extracting all neurons from ROIs, Percentage TH treads are being analysed on the ROI now"


// Analyze percentage of TH Threads on ROIs
print "Start process of pixel classifier. Detect all TH threads in ROI" 
selectAnnotations();
addPixelClassifierMeasurements("TH threads", "TH threads")
saveAnnotationMeasurements('/C:/path/to/output/folder/TH_threads/')
 
 print "Done with whole script"
