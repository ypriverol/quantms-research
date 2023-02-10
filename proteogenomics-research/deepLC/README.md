# DeepLC_Report

In order to use DeepLC to predict the retention time of non-canonical peptides, we did a prediction of canonical peptide to see how the model works. All canonical peptides are divided into ten parts, nine for model calibration and one for prediction. I used two calibration model methods. At the same time, the direct prediction method of uncalibrated model is also used. Here are the results of the predicted canonical peptide.

### pygam_calibration

```python
deeplc --file_pred can1.csv --file_cal can9.csv --plot_predictions
```

![can1_deeplc_predictions](https://user-images.githubusercontent.com/83333440/218061182-2f28b887-e4cc-499e-b261-c118ad8b2902.png)

### legacy_calibration
```python
deeplc --file_pred can1.csv --file_cal can9.csv --plot_predictions --legacy_calibration
```

![can1_deeplc_predictions](https://user-images.githubusercontent.com/83333440/218061300-4d052787-25e4-41d0-add2-fe4faf716cfa.png)
### no calibration
```python
deeplc --file_pred can1.csv --plot_predictions
```

![image](https://user-images.githubusercontent.com/83333440/218063568-0b5c5838-3876-4e46-96f0-0f220c8c40f2.png)

At the same time, I used canonical peptides to calibrate the model to predict non-canonical peptides, and the results are as follows:

### pygam_calibration
```python
 deeplc --file_pred non.csv --file_cal can.csv --plot_predictions
```

![non_to_rt_deeplc_predictions](https://user-images.githubusercontent.com/83333440/218064924-9fceb3d2-eab5-4dab-a3e4-44b34b622c15.png)

### legacy_calibration

```python 
deeplc --file_pred non.csv --file_cal can.csv --plot_predictions --legacy_calibration
```

![non_deeplc_predictions](https://user-images.githubusercontent.com/83333440/218066151-bb1f5f2e-0c35-4dcf-9e41-2671cadb6512.png)

### no calibration

```python 
deeplc --file_pred non.csv --plot_predictions
```

![image](https://user-images.githubusercontent.com/83333440/218065199-7efd37ae-81bc-4677-9c7e-7d09d92bd804.png)