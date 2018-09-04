<!-- dx-header -->
# stjude_chipseq_parameter_wrapper (DNAnexus Platform App)

Wrapper application for ChIP-seq pipeline

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

<!-- Insert a description of your app here -->

# St. Jude ChIP-seq Parameter Wrapper

This app is a setup step for the ChIP-seq pipeline. It configures the pipeline based on the chosen settings.

The output folder is not necessary, and can be left as '/'.
The output prefix **is** required. 

The final outputs will be found in the following path: OUTPUT_FOLDER/Results/OUTPUT_PREFIX. This path will be created by the pipeline. 

> If the above path already exists, make sure it is **empty**. The pipeline deletes intermediate files in the output path and therefore requires an empty directory.
> If the path exists and is non-empty, then a part of the unique job id will be
> added to the folder name.

-->
