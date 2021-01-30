# DELA Instructions from David Lambright

After downloading, double-click "Dela App Scripts.zip" to unzip it. You should get a folder "Dela App Scripts" that contains the application and various other folders/files. You can move the application to your MACs Application folder if desired. It should run on any MAC with OSX 9 or later. Some (admittedly dated) documentation is available within the program under help in the main menu bar. There are also some help tags that can be viewed by hovering the mouse over buttons, etc.

Probably the easiest way to get going would be to try to reproduce the Integrated Michaelis-Menten Analysis for the WT and mutant data you sent me. The steps would be something like:

1. Open the application - should get a new "Untitled" document window. See Document window under the help contents for explanation. You can save the document at any point under the File menu in the toolbar.

2. Import your data - this will likely require conversion into a format that DELA can recognize. Basically, DELA stores data as columns of X, Y, and Errors. It can import data from files or the clipboard as some combination X, Y and/or Errors. In addition to columns of numbers, DELA can import data and/or columns labels separated by a tab. It can also import row labels terminated by a tab. It doesn't know what to do with multiple row labels or formats were one data set is appended after another. That said, it can import data organized as repeating columns in a spread sheet (e.g. X Y X Y etc, or X Y E, X Y E, etc or X Y1 Y2 Y3 etc).

One way to import data would be to open the files in Excel, select/copy the X and Y columns for a given well and paste into DELA. Another way would be to drag/drop files with the X and Y columns for a given well.

After importing the data, a new sheet containing the data will appear in the document contents. The sheet may contain one or more data sets.

3. Plot the data - select the sheet(s) containing the data of interest and plot the data with New Plot under the Plot menu. Hopefully, it looks like your data. Objects in the plots including data points, axes, text, etc are clickable. You can also right click any where on the plot view to get a menu with various options.

4. Select a model - choose Select under the Model menu and select Integrated Michaelis Menten from the available models in the Model popup menu.

5. Setup the initial parameter values for the model - click the Parameters tab above the table and click the Estimate button. This will provide reasonable estimates for the Initial and Final Conversion Factors. You'll need to enter the Initial [Substrate] and Toal [Enzyme] concentrations and make reasonable guesses for kcat and Km. Concentrations units can be anything (e.g. pM, nM, uM, mM, etc) but should be consistent (.e.g all uM or all nM, etc). You can also enter a value for the intrinsic rate (kintr) if known. Assuming your guesses are close enough, you should see model curves that look qualitatively like the data but differ quantitatively. If not, play around with the kcat and Km values. You can use the data popup menu to select the data set of interest if you have more than one data set in the plot. Also, if the Apply to all button is checked, any edits in the table will be applied to all data sets in the plot.

6. Fit the data - choose Fit from the tab above the table and click on either the Marquardt or Simplex buttons. These are the two most common algorithms for fitting models to data. Usually, SImplex is more powerful but sometimes Marquardt works better. I usually alternate between them until the model curves and parameter values stop changing. Probably you will want to fit four parameters: the initial and final conversion factors, kcat, and Km. For the higher concentration data that slopes downward at longer times, you can improve that fit by including a small negative Baseline Slope. There's a Fit column in the table with check boxes determine whether a given parameter is treated as adjustable or as a constant.

Hopefully, thats enough to get started. If you run to problems or have questions, let me know.


## Specifics from Tina

- Import individual files containing three columns: condition, Time, fluorescence
- On import screen, select "has row labels", and set the Data Type as "Fluor Kinetics"
- Tina organizes the datasets in a sheet by copying the mutant descriptor (WT or F28V, etc) before the condition, e.g. the dataset is named "Data4: WT 20201117-WT-A1-5-1-20" (you can't get rid of the Data#). Then you can reorder them by mutant for easy browsing by dragging the datasets
- select a dataset, plot it, then go to model -> IMM
- uncheck all the global fit checkboxes, all of Tina’s fits were done on individual curves
- Only check the following boxes for fitting:
    - initial conversion
    - final conversion
    - baseline slope
    - kcat
    - Km
- If you have bad data, you can try to constrain min and max values in the fit menu, like for the Km
- Initial vals:
   - Put in initial values for kcat 5, Km 0.2
    - Total enzyme 0.001
    - Initial Substrate calculated from the script (file called "calculated_product_conc.txt")
- Running the fit
    - Select Estimate, confirm the global fits are off
    - then fit, Start with simplex, then marquardt. simplex won't fit the baseline
    - Can mask later data points to get better fits of the part that matters, but try not to mask much. You can look and see if that values don’t change too much
- To export:
    - Go to parameters
    - Par -> Data
    - Then go to Sheet: Parameters for IMM -> Parameters for IMM
    - Manually copy into a text file
    - Remove "Data#:" from each line
- To compute errors and mean parameter values, use the "Postdela analysis" file
    - the input to that is the parameters file made by hand from DELA export


## Additional notes from CHris
- to disable dark mode, type `defaults write com.dgl.Dela NSRequiresAquaSystemAppearance -bool yes`, from https://apple.stackexchange.com/questions/338044/can-i-turn-to-dark-mode-only-for-specific-apps-in-macos-mojave