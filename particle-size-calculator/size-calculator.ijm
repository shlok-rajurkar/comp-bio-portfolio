//GGE particle size macro for ImageJ to replace GGE7A3 on NIHImage.

main();


function main() {
    //print('main');
    initialize();
    getStandards();
    cubicRegression(stdRfVals, stdWeights);
    getLanes();
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

function getLanes() {
    //print('getLanes');
    moreLanes = true;
    while (moreLanes) {
        quantifyLane();
        moreLanes = getBoolean("Analyze more lanes?");
    }
}
var stdWeights;
var stdRfVals;
var slopeAndInterceptArray;
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
var yValsCurrLane;
var bins;
var originXValForBinning;


function calcRfVals(xVals, imgLength) {
    //print('calcRfVals');
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = xVals[i]/imgLength;
    }
    return result;
}

function log10Array(array) {
    //print('log10Array');
    logValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        logValues[i] = log(array[i])/log(10);
    }
    return logValues;
}

function inverseLog10Array(array) {
    //print('inverseLog10Array');
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}

function getRfValsFromLane() {
    //print('getRfValsFromLane');

    yValsCurrLane = getProfile();

    run("Plot Profile");

    run("Remove Overlay");

    setTool("multi-point");

    waitForUser("Mark origin of lane, then press OK.");



    run("Measure");

    originXValForBinning = getResult("X");

    getSelectionCoordinates(originXVal, originYVal);

    waitForUser("Mark the peaks with the multipoint tool, then press OK. \nOnly mark them from left to right, in increasing weight/diameter.");
    

    getSelectionCoordinates(xpoints, ypoints);

    nPoints = xpoints.length;

    xVals = newArray(nPoints);

    for (i = 0; i < nPoints; i++) {
        xVals[i] = xpoints[i];
    }

    xVals = Array.slice(xVals, 1, xVals.length);

    getDimensions(plotWidth, plotHeight, c, z, t);

    RfVals = calcRfVals(xVals, plotWidth-originXVal[0]);

    return RfVals;    
}

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

function calcStdCurveLinear(stdRfValues, log10StdWeights) {
    //print('calcStdCurveLinear');
    slopeAndInterceptArray = OLSRegression(stdRfValues, log10StdWeights);  
    return slopeAndInterceptArray;  
}

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



function initialize() {
    //print('initialize');

    imageWidth = 0
    imageHeight = 0
    xnum = ""

    
    run("Set Measurements...", "invert redirect=None decimal=1");
    run("Gel Analyzer Options...", "vertical=1 horizontal=1 label");

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

function cubicRegression(stdRfValues, stdWeights) {
    //print('cubicRegression');
    log10StdWeights = log10Array(stdWeights);
    cubicCoeffArray = calcStdCurveCubic(stdRfVals, log10StdWeights);
    //Array.print(cubicCoeffArray);
    return cubicCoeffArray;
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

function inverseLog10Array(array) {
    //print('inverseLog10Array');
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}


function calcPxFromBins() {
    // Need to fix this, probably where the bin problem is. Bin widths seem to match up with desired widths, but maybe the lateral placement is off...
    //print('calcPxFromBins');
    bins = newArray(
    375, 339, 321, 315, 309, 303, 297, 291, 285, 272, 265, 256, 247, 242, 233, 220
    );
    //bins = newArray(100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20);
    originIndex = originXValForBinning/0.005;
    print(originIndex);
    laneLength = yValsCurrLane.length-originIndex;
    binCount = bins.length; 
    everyPixel = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyPixel[i] = i;
    }
    everyRfValue = newArray(laneLength);
    everyLogMW = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyRfValue[i] = everyPixel[i]/laneLength;
    }

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



function quantBins() {
    //print('quantBins');
    binPxValues = calcPxFromBins();
    baselineY = getBackgroundConc();
    Array.getStatistics(yValsCurrLane, min, max, mean, stdDev);
    print(min, max, mean, stdDev);
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

function sumSingleBin(yValArray, binLowerBound, binUpperBound) {
    //print('sumSingleBin');
    slicedArray = Array.slice(yValArray, binLowerBound, binUpperBound);
    //Array.print(slicedArray);
    Array.getStatistics(slicedArray, min, max, mean, stdDev);
    //print(min, max, mean, stdDev);
    print(slicedArray.length);
    return slicedArray.length * mean;
}



function getBackgroundConc() {
    //print('getBackgroundConc');
    // Can write this in later, should focus on the actual integration right now
    // if(selectingBaseline) {
    //     addHorizontalLineToCursor();
    // }
    // function addHorizontalLineToCursor() {
    //     while(true) {
    //         getCursorLoc(cursorX, cursorY, cursorZ, cursorModifiers);
    //         Overlay.remove;
    //         Overlay.drawLine(0, cursorY, laneLength, cursorY);
    //         Overlay.add;} 
    setTool("multi point");
    waitForUser("Select point that reflects baseline y-value of the LDL range, then press OK.");
    run("Measure");
    baselineY = getResult("Y");

    print(baselineY);
    return baselineY
    }

