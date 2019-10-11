Hierarchical image segmentation provides region-oriented scale-spaces:
sets of image segmentations at different detail levels in which the
segmentations at finer levels are nested with respect to those at
coarser levels. This is the code for Hierarchical Graph Based Segmentation from a non Increasing Edge Observation Attribute. [Article](https://link.springer.com/chapter/10.1007/978-3-030-14085-4_14). 

### Authors






<table style="width:100%; border-collapse: collapse; border: none;">
  <tr>
<td>Edward Cayllahua Cahuina</td>
<td></td>
<td> Université Paris-Est, LIGM, CNRS - ENPC - ESIEE Paris -UPEM </td>
</tr>
<tr>   
<td>Jean Cousty</td>
<td></td>
<td> Université Paris-Est, LIGM, CNRS - ENPC - ESIEE Paris -UPEM </td>
</tr>
<tr>
<td>Silvio Guimarães</td>
<td></td>
<td> Pontifical Catholic University of Minas Gerais, Computer Science Department, VIPLAB</td>
</tr>
<tr>
<td>Yukiko Kenmochi</td>
<td></td>
<td> Université Paris-Est, LIGM, CNRS - ENPC - ESIEE Paris -UPEM </td>
</tr>
<tr>
<td>Guillermo Cámara-Chávez</td>
<td></td>
<td> Universidade Federal de Ouro Preto, Computer Science Department</td>
</tr>
<tr>
<td>Arnaldo de Albuquerque Araújo</td>
<td></td>
<td> Universidade Federal de Minas Gerais, Computer Science Department</td>
</tr>
</table>

### Results
## Video
(TO BE INCLUDED)
## Resulting hierarchy and segmentations
<table style="width:100%">
  <tr>
    <th>Image</th>
    <th>Saliency map for original Min-rule</th>
    <th>Saliency map for proposed Upper P-rank</th>
  </tr>
  <tr>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/2010_000666.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/MinSM.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/UpperPrankSM.png" width="255" height="166"></td>
  </tr>
  <tr>
    <td rowspan="2"></td>
    <th>A segmentation obtained from Min-rule hierarchy</th>
    <th>A segmentation obtained from Upper P-rank hierarchy</th>
  </tr>
  <tr>
    
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/Min_rule.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/upperPrank.png" width="255" height="166"></td>
  </tr>
</table>
Images used for tests come from PASCAL VOC 2010 and VOC 2012 dataset. Files can be downloaded from [LINK](https://drive.google.com/drive/folders/1zYyttKRxeCmk5235wKfV3hfMmeoRSdK2)
### Get the code
The code is available at [CODE](https://github.com/cayllahe/hgbcode).

### Build/Install
This code was compiled and executed in Linux. To compile use:
``` 
./compileHGB.sh
```
### Usage 
The program takes as input an image file in the format XXXX.ppm and produces as output a saliency map XXXX.pgm, which is the visualization of the hierarchy after performing the hierarchical graph based segmentation from a non increasing edge observation attribute of the input image. 

To execute the program: 

```
./hgbSegmentation.sh INPUT_IMAGE.ppm OUTPUT_SM.pgm OPTION TH_PARAMETER RANK_PARAMETER
```
Where:  

**OPTION**:

```
1: Use min-rule: select minimum value on positive observation intervals
2: Use max-rule: select the last upper bound on negative intervals

3: On positive intervals, apply length threshold and min-rule
4: On negative intervals, apply length threshold and max-rule

5: On positive intervals, apply rank filter and  min-rule
6: On negative intervals, apply rank filter and  max-rule

7: On positive intervals, apply area and  min-rule
8: On negative intervals, apply area and  max-rule

9: On negative intervals, apply area, then ranking and max-rule

10: On positive intervals, apply volume filter and min-rule
11: On negative intervals, apply volume filter and max-rule

12: On positive intervals, apply depth filter and min-rule
13: On negative intervals, apply depth filter and max-rule
```

**TH_PARAMETER**: Refers to the  alpha parameter (threshold) 

**RANK_PARAMETER**: Refers to the ranking filtering parameter

**Example**: 
```
./hgbSegmentation.sh Images/3063.ppm /tmp/salida.pgm 11 100 0.1 
```

### Citation 
Please cite as

```
@InProceedings{Cayllahuaetal2019,
author="Cayllahua-Cahuina, Edward
and Cousty, Jean
and Guimar{\~a}es, Silvio
and Kenmochi, Yukiko
and C{\'a}mara-Ch{\'a}vez, Guillermo
and de Albuquerque Ara{\'u}jo, Arnaldo",
title="A Study of Observation Scales Based on Felzenswalb-Huttenlocher Dissimilarity Measure for Hierarchical Segmentation",
booktitle="Discrete Geometry for Computer Imagery",
year="2019",
publisher="Springer International Publishing",
pages="167--179",
isbn="978-3-030-14085-4"
}
```


