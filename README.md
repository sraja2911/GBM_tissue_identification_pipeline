# I. GBM_tissue_identification_pipeline
This pipeline has following patho-image processing methodologies and steps. These steps are implemented using Scikit, HistomicsTK libraries in Python. TCGA-GBM images (around 2000 WSI slides) will be used and DSA platform is used to download them. All the wsi slides will be converted to jpg format, scaled down to the factor of 32. 

1. Image Filtering:
   Filters to Grayscale->RGB to Grayscale
   Complement the grayscaled images
   
2. Adaptive histogram equalization thresholding will be applied on the images - adaptive histogram equalization applies transformations to local regions in an image. As a result, adaptive equalization allows contrast to be enhanced to different extents in different regions based on the regions' intensity histograms. For more information about adaptive equalization, please see https://en.wikipedia.org/wiki/Adaptive_histogram_equalization

3. K-means segmentation on OTSU threasholded images with 3000 segments 
The scikit-image library contains functionality that allows for image segmentation using k-means clustering based on location and color. This allows regions of similarly colored pixels to be grouped together. These regions are colored based on the average color of the pixels in the individual regions. This could potentially be used to filter regions based on their colors, where we could filter on pink shades for eosin-stained tissue and purple shades for hematoxylin-stained tissue.

4. Various color filtering will be applied: Red, Green, Blue, Gray, HSV (Hue-Saturation-Value)

5. Image morphology: Erosin, dilation and Adaptive opening (erosion followed by dilation)
Information about image morphology can be found at https://en.wikipedia.org/wiki/Mathematical_morphology. The primary morphology operators are erosion, dilation, opening, and closing. With erosion, pixels along the edges of an object are removed. With dilation, pixels along the edges of an object are added. Opening is erosion followed by dilation. Closing is dilation followed by erosion. With morphology operators, a structuring element (such as a square, circle, cross, etc) is passed along the edges of the objects to perform the operations. Morphology operators are typically performed on binary and grayscale images.

6. Entropy filtering

7. Tiles or ROI identified 


# II. Necrotic and MicroVascular proliferation (MVP) identification  pipeline

All the above identified tissues/tiles will be used for

1. An automatic cell segmentation, with cell-profile count with Decision tree will be implemented for the Necrosis ("cellular Pseudopalisades") identification
2. Spatial image filtering, statistical techniques on MVP identification  
DSA_viewer plug-in will be used for ground truth annotations on these images. Histoclassifier (TensorFlow based Layer multi classifier, will be customized to suit).

Reference:
1. http://cv-tricks.com/tensorflow-tutorial/training-convolutional-neural-network-for-image-classification/ - "A TensorFlow binary classifier, modified to multi-classifier"

2. https://github.com/CODAIT/deep-histopath/blob/master/docs/wsi-preprocessing-in-python/index.md - Python code for WSI preprocessing - "Customized version for GBM images and training path and HistomicTK utilites"

3.Automated discrimination of lower and higher grade gliomas based on histopathological image analysis, Hojjat Seyed Mousavi, Vishal Monga, Ganesh Rao, Arvind U. K. Rao - Jan'15

4. 'Pseudopalisading' Necrosis in Glioblastoma: A familiar Morphologic Feature That Links Vascular Pathology, Hypoxia and Angiogenesis, Yuan Rong, Donald L.Durden, Erwin G.Van Meir and Danier J.Brat - June'06

5. The Tumor Microenvironment Strongly Impacts Master Transcriptional Regulators and Gene Expression Class of Glioblastoma, Lee A.D. Cooper, David A. Gutman, Candace Chisolm, Christina Appin, Jun Kong, Yuan Rong,Tahsin Kurc, Erwin G. Van Meir,Joel H. Saltz,Carlos S. Moreno and Daniel J. Brat - May'12

6. The 2007 WHO Classification of Tumours of the Central Nervous System,David N. Louis · Hiroko Ohgaki · Otmar D. Wiestler, Webster K. Cavenee · Peter C. Burger · Anne Jouvet, Bernd W. Scheithauer · Paul Kleihues, May'07
