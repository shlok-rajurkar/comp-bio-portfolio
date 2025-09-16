//GGE particle size macro for ImageJ to replace GGE7A3 on NIHImage.

main();

//Global vars, will transition away from global vars once working properly.

var stdWeights;
var stdRfVals;
var imageWidth;
var imageHeight;
var xnum;
var c = 0;
var z = 0;
var t = 0;
var cubicCoeffArray;
var croppedGelWindow;
var originXVal;
var laneLength;
var xValsCurrLaneIndex;
var yValsCurrLane;
var bins;

// Main macro flow
function main() {
    //print('main');
    initialize();
    getStandards();
    cubicRegression(stdRfVals, stdWeights);
    Array.print(cubicCoeffArray);
    getLanes();
}

// Prompts user to crop gel
function initialize() {
    //print('initialize');

    imageWidth = 0
    imageHeight = 0
    xnum = ""

    
    run("Set Measurements...", "invert redirect=None decimal=3");
    run("Gel Analyzer Options...", "vertical=1 horizontal=1 label");
    run("Clear Results");

    xnum = getString("enter X number", "");
    
    setTool("Rectangle");

    getDimensions(imageWidth, imageHeight, c, z, t);
    makeRectangle(0.3*imageWidth, 0.18*imageHeight, 0.36*imageWidth, 0.53*imageHeight);
    waitForUser("Adjust rectangle to bound lanes and leave minimal space between top of image and top of gel. \nPressing OK will crop gel and rotate 90 degrees left.");

    run("Crop");
    wait(100);
    run("Rotate 90 Degrees Left");
    croppedGelWindow = File.name;
    print(xnum + " Analysis:");
}

function getStandards() {
    //print('getStandards');
    stdWeights = newArray(0);
    stdRfVals = newArray(0);
    moreStandards = true;
    while (moreStandards) {
        stdWeightsTempAndStdRfValsTemp = setStandards();
        arrayLength = stdWeightsTempAndStdRfValsTemp.length;
        stdWeightsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, 0, (arrayLength/2));
        stdRfValsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, (arrayLength/2), arrayLength);
        stdWeights = Array.concat(stdWeights, stdWeightsToAdd);
        stdRfVals = Array.concat(stdRfVals, stdRfValsToAdd);
        stdWeights = Array.sort(stdWeights);
        stdWeights = Array.reverse(stdWeights);
        stdRfVals = Array.sort(stdRfVals);
        moreStandards = getBoolean("Add more standard lanes?");
    }
}

function setStandards() {
    //print('setStandards');
    selectWindow(croppedGelWindow);
    numberOfStandards = getNumber("enter number of standards", 5);
    stdWeightsTemp = newArray(numberOfStandards);
    stdRfValsTemp = newArray(numberOfStandards);
    for(i = 0; i < numberOfStandards; i++){
        stdIndexDisplay = i + 1;
        stdWeightsTemp[i] = getNumber("enter weight of standard " + stdIndexDisplay, 0);
    }
    print("Std weights:");

    setTool("Rectangle");

    getDimensions(standardLaneWidth, standardLaneHeight, c, z, t);

    makeRectangle(0, (standardLaneHeight/2)-10 , standardLaneWidth, 20);
    waitForUser("Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.");
    
    run("Select First Lane");

    stdRfValsTemp = getRfValsFromLane();
    stdWeightsTempAndStdRfValsTemp = Array.concat(stdWeightsTemp, stdRfValsTemp);
    return stdWeightsTempAndStdRfValsTemp;
}

function getRfValsFromLane() {
    //print('getRfValsFromLane');

    run("Plot Profile");

    Plot.getValues(xValsCurrLane, yValsCurrLane);
    // Don't think subtracting index 0 is necessary here since X vals will always start at 0.
    xValsCurrLaneIndex = divideArrayByStep(xValsCurrLane, xValsCurrLane[1]-xValsCurrLane[0]);

    run("Remove Overlay");

    setTool("multi-point");

    waitForUser("Mark origin of lane, then press OK.");

    run("Clear Results");

    run("Measure");

    originXVal = divideByStep(getResult("X"), xValsCurrLane[1]-xValsCurrLane[0]);

    laneLength = xValsCurrLaneIndex.length - originXVal;

    xValsCurrLaneIndex = Array.slice(xValsCurrLaneIndex, originXVal, xValsCurrLaneIndex.length);

    for (i = 0; i < xValsCurrLaneIndex.length; i ++) {
        xValsCurrLaneIndex[i] = xValsCurrLaneIndex[i] - originXVal;
    }

    yValsCurrLane = Array.slice(yValsCurrLane, originXVal, yValsCurrLane.length);

    run("Clear Results");

    run("Remove Overlay");

    waitForUser("Mark peaks with the multipoint tool, then press OK. \nOnly mark them from left to right, in increasing weight/diameter.");

    run("Measure");

    xValsAndOrigin = divideArrayByStep(getAllResults("X"), xValsCurrLane[1]-xValsCurrLane[0]);

    Array.print(xValsAndOrigin);

    xVals = Array.slice(xValsAndOrigin, 1);

    Array.print(xVals);

    RfVals = calcRfVals(xVals, laneLength);

    return RfVals;    
}

getRfValsFromLaneLineCursor() {
    return None
}

function getLanes() {
    //print('getLanes');
    moreLanes = true;
    while (moreLanes) {
        quantifyLane();
        moreLanes = getBoolean("Analyze more lanes?");
    }
}

function quantifyLane() {
    //print('quantifyLane');
    if (stdRfVals.length == 0 || stdWeights.length == 0){
        exit("Standards not set.");
    }

    selectWindow(croppedGelWindow);

    waitForUser("Drag rectangle to sample lane and press OK.");

    run("Select First Lane");

    laneRfVals = getRfValsFromLane();

    laneMolecularWeightsCalc = newArray(laneRfVals.length);
    for (i = 0; i < laneRfVals.length; i++) {
        x = laneRfVals[i];
        laneMolecularWeightsCalc[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*x + cubicCoeffArray[2]*x*x + cubicCoeffArray[3]*x*x*x;
    }
    displayValueArray = inverseLog10Array(laneMolecularWeightsCalc);
    print("Calculated Diameters:");
    Array.print(displayValueArray);

    quantBins();

}

function quantBins() {
    //print('quantBins');
    binPxValues = calcPxFromBins();
    baselineY = getBackgroundConc();
    Array.getStatistics(yValsCurrLane, min, max, mean, stdDev);
    //print(min, max, mean, stdDev);
    //Array.print(yValsCurrLane);
    for (i = 0; i < yValsCurrLane.length; i++) {
    yValsCurrLane[i] = yValsCurrLane[i] - baselineY;
    }
    binSums = newArray(binPxValues.length-1);
    for (i = 0; i < binPxValues.length-1; i ++) {
        binSums[i] = sumSingleBin(yValsCurrLane, binPxValues[i], binPxValues[i+1]);
    }
    Array.getStatistics(binSums, min, max, mean, stdDev);
    binSumsTotal = mean*binSums.length;
    for (i = 0; i < binSums.length; i++) {
        binSums[i] = (binSums[i]/binSumsTotal)*100;
        print(bins[i] + '-' + bins[i+1] + ": " + binSums[i] + "%");
    }
    
    //Array.print(binSums);

}

function calcPxFromBins() {
    //print('calcPxFromBins');
    bins = newArray(
    375, 339, 321, 315, 309, 303, 297, 291, 285, 272, 265, 256, 247, 242, 233, 220
    );
    binCount = bins.length; 
    everyRfValue = calcRfVals(xValsCurrLaneIndex, laneLength);
    everyLogMW = newArray(xValsCurrLaneIndex.length);
    for (i = 0; i < laneLength; i ++) {
            everyLogMW[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*everyRfValue[i] + cubicCoeffArray[2]*everyRfValue[i]*everyRfValue[i] + cubicCoeffArray[3]*everyRfValue[i]*everyRfValue[i]*everyRfValue[i];
    }
    everyMW = inverseLog10Array(everyLogMW);
    binPxValues = newArray(binCount);
    for (j = 0; j < bins.length; j ++) {
        everyMWCopy = Array.copy(everyMW);
        for (i = 0; i < laneLength; i ++) {
            everyMWCopy[i] = abs(everyMWCopy[i] - bins[j]);   
        }
        Array.getStatistics(everyMWCopy, min, max, mean, stdDev);

        for (h = 0; h < laneLength; h ++) {
            if (abs(everyMWCopy[h] - min)<1e-6) {
                binPxValues[j] = h;
            }
        }
    }
    //Array.print(binPxValues);
    return binPxValues;
}

function sumSingleBin(yValArray, binLowerBound, binUpperBound) {
    //print('sumSingleBin');
    slicedArray = Array.slice(yValArray, binLowerBound, binUpperBound);
    //Array.print(slicedArray);
    Array.getStatistics(slicedArray, min, max, mean, stdDev);
    //print(min, max, mean, stdDev);
    //Array.print(slicedArray);
    return slicedArray.length * mean;
}

function cubicRegression(stdRfValues, stdWeights) {
    //print('cubicRegression');
    log10StdWeights = log10Array(stdWeights);
    cubicCoeffArray = calcStdCurveCubic(stdRfVals, log10StdWeights);
    //Array.print(cubicCoeffArray);
    return cubicCoeffArray;
}

function getBackgroundConc() {
    //print('getBackgroundConc');

    setTool("multi point");
    waitForUser("Select point that reflects baseline y-value of the LDL range, then press OK.");
    run("Clear Results");
    run("Measure");
    baselineY = getResult("Y");

    //print(baselineY);
    return baselineY;
}



// -- Utils -- //

// Rounds number to given decimal step
function roundToStep(number, step) {
    return round(number/step)*step;
}

// roundToStep for array
function roundArrayToStep(array, step) {
    roundedArray = newArray(array.length);
    for (i = 0; i < array.length; i ++) {
        roundedArray[i] = roundToStep(array[i], step);
    }
    return roundedArray;
}

// Returns whole multiples of step in number
function divideByStep(number, step) {
    return round(number/step);
}

// divideByStep for array
function divideArrayByStep(array, step) {
    dividedArray = newArray(array.length);
    for (i = 0; i < array.length; i ++) {
        dividedArray[i] = divideByStep(array[i], step);
    }
    return dividedArray;
}

// Retrieves all results currently in Results window
function getAllResults(column) {
    array = newArray(nResults);
    results = newArray(array.length);
    for (i = 0; i < results.length; i ++) {
        results[i] = getResult(column, i);
    }
    return results;
}

// Calculates Rf vaules given an array of x values and total lane length
function calcRfVals(xVals, laneLength) {
    //print('calcRfVals');
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = xVals[i]/laneLength;
    }
    return result;
}

// Calculates log base 10 for every value of an array
function log10Array(array) {
    //print('log10Array');
    logValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        logValues[i] = log(array[i])/log(10);
    }
    return logValues;
}

// Calculates inverse log 10 for every value of an array
function inverseLog10Array(array) {
    //print('inverseLog10Array');
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}

// Performs OLS regression for arrays of x values and y values
function OLSRegression(xValsForRegression, yValsForRegression) {
    //print('OLSRegression');
    xLength = xValsForRegression.length;
    yLength = yValsForRegression.length;

    xTotal = 0;
    yTotal = 0;

    for(i = 0; i < xLength; i++) {
        xTotal += xValsForRegression[i];
        yTotal += yValsForRegression[i];
    }
    xMean = xTotal/xLength;
    yMean = yTotal/yLength;

    numerator = 0;
    denominator = 0;

    for(i = 0; i < xLength; i++) {
        dx = xValsForRegression[i] - xMean;
        dy = yValsForRegression[i] - yMean;
        numerator += dx * dy;
        denominator += dx * dx;
    }

    regressionSlope = numerator/denominator;

    regressionIntercept = yMean - regressionSlope * xMean;

    return newArray(regressionSlope, regressionIntercept);
}

// Calculates cubic regression coefficients given x values and y values
function calcStdCurveCubic(xVals, yVals) {
    //print('calcStdCurveCubic');
    n = xVals.length;

    sumX = newArray(7);
    for (i = 0; i <= 6; i++) sumX[i] = 0;
    sumXY = newArray(4);
    for (i = 0; i <= 3; i++) sumXY[i] = 0;

    for (i = 0; i < n; i++) {
        x = xVals[i];
        y = yVals[i];
        powX = newArray(7);
        powX[0] = 1;
        for (j = 1; j <= 6; j++) powX[j] = powX[j - 1] * x;

        for (j = 0; j <= 6; j++) sumX[j] += powX[j];
        for (j = 0; j <= 3; j++) sumXY[j] += powX[j] * y;
    }

    // Flattened 4x4 matrix A
    A = newArray(16);
    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++) {
            A[i * 4 + j] = sumX[i + j];
        }
    }

    B = sumXY;

    coeffs = gaussJordanFlat4x4(A, B);
    return coeffs; // [a, b, c, d]
}

// Gauss Jordan reduction of 4x4 matrix
function gaussJordanFlat4x4(A, B) {
    //print('gaussJordanFlat4x4');
    n = 4;

    for (i = 0; i < n; i++) {
        // Find non-zero pivot
        if (A[i * 4 + i] == 0) {
            for (j = i + 1; j < n; j++) {
                if (A[j * 4 + i] != 0) {
                    // Swap rows in A
                    for (k = 0; k < n; k++) {
                        temp = A[i * 4 + k];
                        A[i * 4 + k] = A[j * 4 + k];
                        A[j * 4 + k] = temp;
                    }
                    // Swap B
                    tmpB = B[i];
                    B[i] = B[j];
                    B[j] = tmpB;
                    break;
                }
            }
        }

        // Normalize row
        factor = A[i * 4 + i];
        for (j = 0; j < n; j++) A[i * 4 + j] /= factor;
        B[i] /= factor;

        // Eliminate other rows
        for (j = 0; j < n; j++) {
            if (j != i) {
                factor = A[j * 4 + i];
                for (k = 0; k < n; k++) {
                    A[j * 4 + k] -= factor * A[i * 4 + k];
                }
                B[j] -= factor * B[i];
            }
        }
    }

    return B;
}

// Connects cursor to horizontal line
function addHorizontalLineToCursor() {
    getDimensions(w, h, c, z, f);
    while(true) {
        getCursorLoc(x, y, z, m);
        Overlay.drawLine(0, y, w, y);
        Overlay.show();
        wait(5);
        Overlay.remove();
}
}

// Problem now is how to have the wait for user function coexist with the vert/horiz line code blocks. 
// I don't think I can run the wait for user without exiting the code block. Might be cooked, but glad
// I could implement the function in the first place. I don't think I can maintain the ability to undo 
// either if I make the function quit after a certain number of points are selected. I would probably 
// need some way for the user to signal that they have selected all points while remaining in the loop.
// Some var that can be checked during the loop... One funny way would be detecting whether the window
// is open, having the user close the window, and then having it automatically reopen for the next
// selection. Unsure what other options there are. Need something the user can interact with while they
// are selecting ROIs. This might be the best option unfortunately. Will test it out to see if the small
// change in interactivity has negative effects that are worth the improvement in UI and cursor visibility.
// Suspect that it might be, given that marking the peaks precisely is probably one of the most important 
// factors in a successful analysis.

// Connects cursor to vertical line
function addVerticalLineToCursor() {
    getDimensions(w, h, c, z, f);
    while(true) {
        getCursorLoc(x, y, z, m);
        Overlay.drawLine(x, 0, x, h);
        Overlay.show();
        wait(5);
        Overlay.remove();
    }
}