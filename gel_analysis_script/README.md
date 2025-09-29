# PAGE Particle Size Calculator
This is a script usable in ImageJ to analyze the particle size of bands in polyacrylamide gel electrophoresis (PAGE) images/scans. It supports multiple standard lanes and uses quartic regression. By default, it will calculate the amount of protein found in different LDL bins (determined by diameter). A recent version of ImageJ is required to run this macro and can be found <a href='https://imagej.net/ij/download.html'>here</a>.<br />


Instructions for installation:
1. Open <a href='/particle-size-calculator/size_calculator.ijm' target='_blank'>the size calculator file</a> and click on the download button shown below. Then, return to this page by clicking the back button on your browser. <br />
<img title='download instructions' alt='download instructions' src='/images/download_button_directions.png'>
2. Open ImageJ. <br />
3. In ImageJ, go to Plugins>Macros and select Install. <br />
4. In the File Explorer Window, navigate to the location of size_calculator.ijm (probably in Downloads) and select it. Press Open. <br />
5. The macro can now be initiated from Plugins>Macros>size_calculator. <br />
<br />
<br />

Instructions for use: <br>
**IMPORTANT: Read all instructions carefully. If you receive an error, take note of the error output and restart from step 2.**
1. Open ImageJ.
2. Open gel file via File>Open and navigate to gel image.
3. Initiate macro via Plugins>Macros>size_calculator. If the macro cannot be found, follow installation instructions above.
4. Bound gel with rectangle, leaving a small border so as to not cut off any parts of the gel. You will have a chance later to more precisely mark the origin of each lane. Press OK.
5. For each standard lane:
    1. Enter number of standards.
    1. Enter standard weights.
    1. Position rectangle to capture entire length of lane. It can be quite thin. Press OK. Press 'Yes' if asked whether the lanes are horizontal.
    1. Read displayed instructions for marking origin and standard peaks and press OK. **It is important that the first point you mark is the origin of the lane.**
    1. Select the origin. The origin will often present as a sharp peak signifying the edge of the gel. If you make a mistake, use Alt+Left Click to remove a point.
    1. Mark each standard peak and press the space bar once finished. Make sure the number of standard peaks you mark matches the number of standards you inputted for that lane.
6. For each sample lane:
    1. Position rectangle to capture entire length of lane. It can be quite thin. Press OK.
    1. Read displayed instructions for marking origin and unknowns and press OK.
    1. Select the origin.
    1. Select each unknown and press the space bar once finished.
    1. Read displayed instructions for marking baseline Y value and press OK.
    1. Select baseline Y value and press the space bar once finished.
    1. Continue analyzing lanes if desired.
7. The output in the log window will display the calculated MW/diameters of your selected peaks. By default, the macro will also try to calculate the amount of LDL in each LDL bin. Ignore this data if it is not relevant.

