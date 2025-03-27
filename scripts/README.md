<h1>Running the pipeline</h1>

<p>
 The first step in the workflow running is to obtain variants from WGS and RNA-seq. Both pipelines for these steps can be found in the /BASH scripts. The resluts are VCF files. <br>
 Then after obtaining the read count matrix  R/ script can be applied to define the interesting cases where at least one of the family memebrs has log2FC > 1 and p-adj value < 0.05.<br>
 We achieved our read Count matrix using FeatureCount software, the matrix can be found in the table/ folder. 
</p> 



![Alt text](Slide4.jpg)


