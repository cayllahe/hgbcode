Hierarchical image segmentation provides region-oriented scale-spaces:
sets of image segmentations at different detail levels in which the
segmentations at finer levels are nested with respect to those at
coarser levels. This is the code gfor Hierarchical Graph Based Segmentation from a non Increasing Edge Observation Attribute. [link](https://link.springer.com/chapter/10.1007/978-3-030-14085-4_14). 

### Images and Results

![oy Image](https://cayllahe.github.io/hgbcode/assets/Figures/2010_000666.png?v=4&s=200)
![another image]({{ site.url }}/assets/2010_000666.png)
![image2]({{ site.url }}/docs/assets/Figures/MinSM.png)
![image3]({{ site.url }}/docs/assets/Figures/UpperPrankSM.png)
![image4]({{ site.url }}/docs/assets/Figures/upperPrank.png)
![image5]({{ site.url }}/docs/assets/Figures/Min_rule.png)

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


