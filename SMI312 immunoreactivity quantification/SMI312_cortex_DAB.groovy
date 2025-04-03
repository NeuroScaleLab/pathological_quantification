// This script quantifies SMI312 immunoreactivity (% area), which targets phosphorylated neurofilament medium and heavy chain

// Script made by Irene Frigerio (2021)
// Used in ex. Frigerio et al. 2023, Translational Neurodegeneration PMID: 36658627
// QuPath version 0.3.2
// Microscope: whole-slide scanner Olympus VS200 Evident, 20x objective
// Antibody clone SMI312 (BioLegend)
// Please do not forget to adjust the output folder


// Standard settings, set deconvolutions for DAB and hematoxylin
setImageType('BRIGHTFIELD_OTHER');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759 ", "Background" : " 255 255 255 "}');


//Run pixel classifier in the ROIs
selectAnnotations();
addPixelClassifierMeasurements("SMI312_percentage_area", "SMI312_percentage_area")
saveAnnotationMeasurements('/C:/path/to/output/folder/')
