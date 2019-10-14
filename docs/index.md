Hierarchical image segmentation provides region-oriented scale-spaces:
sets of image segmentations at different detail levels in which the
segmentations at finer levels are nested with respect to those at
coarser levels. This is the code for Hierarchical Graph Based Segmentation from a non Increasing Edge Observation Attribute <a href="https://link.springer.com/chapter/10.1007/978-3-030-14085-4_14" target="_blank"> [Cayllahua-et-al, 2019] </a>.  

### Authors
<style>
.tablelines table{
        width: 120%;
        border: none!important;
        border-collapse: collapse;       
        border-spacing: 0;        
        }
.tablelines td{
    border: none!important;
    border-collapse: collapse;
}
</style>
<table class="tablelines">
  <tr >
    <td width="35%">Edward Cayllahua-Cahuina [1,2]</td>
    <td width="35%">Jean Cousty [1]</td>
    <td width="35%">Silvio Guimarães [3]</td>    
  </tr>
  <tr>
    <td>Yukiko Kenmochi [1]</td>
    <td>Guillermo Cámara-Chávez[4]</td>
    <td>Arnaldo de Albuquerque Araújo[2]</td>    
  </tr>
    <tr>
    <td colspan="3"> 
        [1] LIGM - ESIEE Paris - Université Paris-Est <br>
        [2] NPDI - DCC - Universidade Federal de Minas Gerais<br>
        [3] VIPLAB - Pontifícia Universidade Católica De Minas Gerais <br>
        [4] DCC - Universidade Federal de Ouro Preto
    </td>    
    </tr>
    <tr> 
    <td colspan="3">    
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/esiee.png" width="90" height="45">
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/ufmg.png" width="90" height="40">
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/ligm.png" width="80" height="70">
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/npdi.png" width="80" height="100">
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/puc.png" width="85" height="70">
        <img src="https://cayllahe.github.io/hgbcode/assets/logos/ufop.png" width="70" height="90">
        </td>   
    </tr>
</table>



### Results
#### Video
<video width="320" height="200" controls preload> 
    <source src="https://github.com/cayllahe/hgbcode/blob/master/docs/assets/video/videoHGB.mp4">     
</video>


#### Resulting hierarchy and segmentations
<table style="width:100%">
  <tr>
    <th>Image</th>
    <th>Result of the original <a href="https://www.degruyter.com/downloadpdf/j/mathm.2017.2.issue-1/mathm-2017-0004/mathm-2017-0004.pdf" target="_blank">HGB method</a> (saliency map)</th>
    <th>Result of the HGB method with the newly proposed  upper P-rank selection strategy</th>
  </tr>
  <tr>
    <td rowspan="3"><img src="https://cayllahe.github.io/hgbcode/assets/Figures/2010_000666.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/MinSM.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/UpperPrankSM.png" width="255" height="166"></td>
  </tr>
  <tr>    
    <th>A segmentation obtained from Min-rule hierarchy</th>
    <th>A segmentation obtained from Upper P-rank hierarchy</th>
  </tr>
  <tr>    
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/Min_rule.png" width="255" height="166"></td>
    <td><img src="https://cayllahe.github.io/hgbcode/assets/Figures/upperPrank.png" width="255" height="166"></td>
  </tr>
</table>
Images used for tests come from PASCAL VOC 2010 and VOC 2012 dataset. Click
[here](https://github.com/cayllahe/hgbcode/tree/master/docs/assets/files) to see more hierarchical image segmentation results 
### Get the code
Click <a href="https://github.com/cayllahe/hgbcode" target="_blank">here</a> to get the code.
### Build/Install
This code was compiled and executed in Linux. To compile use:
``` 
./compileHGB.sh
```
### Usage 
The program takes as input an image file in the format XXXX.ppm and produces as output a saliency map XXXX.pgm, which is the visualization of the hierarchy after performing the hierarchical graph based segmentation from a non increasing edge observation attribute of the input image. 

To execute the program: 

```
./hgbSegmentation.sh INPUT_IMAGE.ppm OUTPUT_SM.pgm OPTION PARAMETER AREASIMP
```
Where:  

**OPTION**:

```
1: Use Min: select minimum value on positive observation intervals
2: Use Max: select the last upper bound on negative intervals

3: Lower-length: On positive intervals, apply length threshold and min-rule
4: Upper-length: apply length threshold and max-rule

5: Lower-area: On positive intervals, apply area and  min-rule
6: Upper-Narea: On negative intervals, apply area and  max-rule

7: Lower-depth: On positive intervals, apply depth filter and min-rule
8: Upper-Ndepth: On negative intervals, apply depth filter and max-rule

9: Lower p-rank: On positive intervals, apply rank filter and  min-rule
10: Upper p-rank: On negative intervals, apply rank filter and  max-rule
```

**PARAMETER**: Refers to the  alpha or ranking parameter (threshold)

**AREASIMP**: Refers to the area simplification parameter. If no value is passed, then the value 0.004 is assumed. 

**Example**: 
```
./hgbSegmentation.sh Images/3063.ppm /tmp/salida.pgm 10 0.001 0.003 
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
### Related articles to HGB

* Efficient Algorithms for Hierarchical Graph-Based Segmentation Relying on the Felzenszwalb–Huttenlocher Dissimilarity <a href="https://repositorio.ufop.br/bitstream/123456789/11338/1/ARTIGO_EfficientAlgorithmsHierarchical.pdf" target="_blank"> [Cayllahua-et-al, 2019] </a>.<br>

* Hierarchizing graph-based image segmentation algorithms relying on region dissimilarity <a href="https://www.degruyter.com/downloadpdf/j/mathm.2017.2.issue-1/mathm-2017-0004/mathm-2017-0004.pdf" target="_blank"> [Guimaraes-et-al, 2017] </a>.<br>

* Hierarchical segmentations with graphs: quasi-flat zones, minimum spanning trees, and saliency maps <a href="https://hal.archives-ouvertes.fr/hal-01344727v2/document" target="_blank"> [Cousty-et-al, 2018] </a>.
