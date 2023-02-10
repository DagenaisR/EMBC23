# EMBC23
Codes used in Dagenais R., Mitsis G. D., Non-invasive estimation of arterial blood pressure fluctuations using a
peripheral photoplethysmograph inside the MRI scanner. EMBC 2023

Processing workflow is shown in Fig. 1




![processing_workflow](https://user-images.githubusercontent.com/102877412/218179480-c493d845-3924-4ab6-a49b-3ad48f07407b.png)
Figure 1. The PPG recordings and the ABP recordings (present for preScan and postScan only) are inputted one at a time inside preProcessing. The output of preProcessing is fed into preModeling. The resulting file is used to train the model used to predict the MAP fluctuations from the PPG signal only (i.e for rest, breathing, and coldPressor). The final output is a matlab binary file containing the desired model-predicted MAP fluctuations.

The instructions are detailed in the codes.
