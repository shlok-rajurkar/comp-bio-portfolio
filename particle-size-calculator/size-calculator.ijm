//GGE particle size macro for ImageJ to replace GGE7A3 on NIHImage.

main();


function main() {
    initialize();
    getStandards();
    cubicRegression(stdRfVals, stdWeights);
    getLanes();
    }


function getStandards() {
    stdWeights = newArray(0);
    stdRfVals = newArray(0);
    moreStandards = true;
    while (moreStandards) {
        stdWeightsTempAndStdRfValsTemp = setStandards();
        arrayLength = stdWeightsTempAndStdRfValsTemp.length;
        stdWeightsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, 0, (arrayLength/2));
        stdRfValsToAdd = Array.slice(stdWeightsTempAndStdRfValsTemp, (arrayLength/2), arrayLength);
        Array.print(stdWeightsToAdd);
        Array.print(stdRfValsToAdd);
        stdWeights = Array.concat(stdWeights, stdWeightsToAdd);
        stdRfVals = Array.concat(stdRfVals, stdRfValsToAdd);
        stdWeights = Array.sort(stdWeights);
        stdWeights = Array.reverse(stdWeights);
        stdRfVals = Array.sort(stdRfVals);
        Array.print(stdWeights);
        Array.print(stdRfVals);
        moreStandards = getBoolean("Add more standard lanes?");
    }
}

function getLanes() {
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

function calcRfVals(xVals, imgLength) {
    result = newArray(xVals.length);
    for (i = 0; i < xVals.length; i ++) {
        result[i] = xVals[i]/imgLength;
    }
    return result;
}

function log10Array(array) {
    logValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        logValues[i] = log(array[i])/log(10);
    }
    return logValues;
}

function inverseLog10Array(array) {
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}

function getRfValsFromLane() {
    run("Plot Lanes");

    run("Remove Overlay");

    setTool("multi-point");

    waitForUser("Mark origin of lane, then press OK.");

    getSelectionCoordinates(originXVal, originYVal);

    waitForUser("Mark the peaks with the multipoint tool, then press OK. \nOnly mark them from left to right, in increasing weight/diameter.");

    getSelectionCoordinates(xpoints, ypoints);

    nPoints = xpoints.length;

    xVals = newArray(nPoints);

    for (i = 0; i < nPoints; i++) {
        xVals[i] = xpoints[i];
    }

    xVals = Array.slice(xVals, 1, xVals.length);

    getDimensions(imageWidth, imageHeight, c, z, t);

    RfVals = calcRfVals(xVals, imageWidth-originXVal[0]);

    return RfVals;    
}

function OLSRegression(xValsForRegression, yValsForRegression) {
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
    slopeAndInterceptArray = OLSRegression(stdRfValues, log10StdWeights);  
    return slopeAndInterceptArray;  
}

function calcStdCurveCubic(xVals, yVals) {
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


    //stdWeights = newArray(0);
    //stdRfVals = newArray(0);
    //slopeAndInterceptArray = newArray(0);
    imageWidth = 0
    imageHeight = 0
    xnum = ""

    run("Set Measurements...", "redirect=None decimal=1");
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
    selectWindow(croppedGelWindow);
    numberOfStandards = getNumber("enter number of standards", 5);
    stdWeightsTemp = newArray(numberOfStandards);
    stdRfValsTemp = newArray(numberOfStandards);
    for(i = 0; i < numberOfStandards; i++){
        stdIndexDisplay = i + 1;
        stdWeightsTemp[i] = getNumber("enter weight of standard " + stdIndexDisplay, 0);
    }
    print("Std weights:");
    Array.print(stdWeightsTemp);

    setTool("Rectangle");

    getDimensions(imageWidth, imageHeight, c, z, t);

    makeRectangle(0, (imageHeight/2)-10 , imageWidth, 20);
    waitForUser("Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.");

    //usePresetRectangle = getBoolean('Use ROI preset?');

    //if (usePresetRectangle) {
    //    makeRectangle(0, (imageHeight/2)-10 , imageWidth, 20);
    //    waitForUser("Adjust rectangle to span lane with standards. \nIt can be quite thin as long as it contains some part of the lane.");
    //}
    //else {
    //    waitForUser("Draw rectangle from left to right of lane with standards. \nIt can be quite thin as long as it contains some part of the lane. Pressing OK will select this lane.");
    //}

    
    run("Select First Lane");

    stdRfValsTemp = getRfValsFromLane();
    Array.print(stdRfValsTemp);
    stdWeightsTempAndStdRfValsTemp = Array.concat(stdWeightsTemp, stdRfValsTemp);
    return stdWeightsTempAndStdRfValsTemp;

    
}

function cubicRegression(stdRfValues, stdWeights) {
    log10StdWeights = log10Array(stdWeights);
    cubicCoeffArray = calcStdCurveCubic(stdRfVals, log10StdWeights);
    Array.print(cubicCoeffArray);
    return cubicCoeffArray;
}

function quantifyLane() {
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
    print("Calculated:");
    Array.print(displayValueArray);

    quantBins() {
        
    }
    originXVal = 0;

}

function inverseLog10Array(array) {
    inverseLogValues = newArray(array.length);
    for (i = 0; i < array.length; i++){
        inverseLogValues[i] = pow(10, array[i]);
    }
    return inverseLogValues;
}
//Test cases for XG4513/14 where stds are 100, 80, 60, 40, 20
imageWidth = 1500;
originXVal = 100;
cubicCoeffArray = newArray(2.0490, 0.0341, -0.8621, -0.1735);
print("-------------------------------------------------------");
function calcPxFromBins() {
    bins = newArray(
    375, 339, 321, 315, 309, 303, 297, 291, 285, 272, 265, 256, 247, 242, 233, 220
    );
    bins = newArray(96, 111, 100, 96, 37, 28);
    laneLength = imageWidth-originXVal;
    binCount = bins.length; 
    everyPixel = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyPixel[i] = i;
    }
    //print('every pixel:');
    //Array.print(everyPixel);
    everyRfValue = newArray(laneLength);
    everyLogMW = newArray(laneLength);
    for (i = 0; i < laneLength; i ++) {
        everyRfValue[i] = everyPixel[i]/laneLength;
    }
    //print('every Rf value:');
    //Array.print(everyRfValue);
    for (i = 0; i < laneLength; i ++) {
            everyLogMW[i] = cubicCoeffArray[0] + cubicCoeffArray[1]*everyRfValue[i] + cubicCoeffArray[2]*everyRfValue[i]*everyRfValue[i] + cubicCoeffArray[3]*everyRfValue[i]*everyRfValue[i]*everyRfValue[i];
    }
    //print('everyMW');
   
    everyMW = inverseLog10Array(everyLogMW);
    //Array.print(everyMW);
    binPxValues = newArray(binCount);
    //Array.print(binPxValues);
    for (j = 0; j < bins.length; j ++) {
        //print("original MWs");
        //Array.print(everyMW);
        everyMWCopy = Array.copy(everyMW);
        //everyMWCopy = everyMW;
        // print("before subtraction:");
        // Array.print(everyMWCopy);
        for (i = 0; i < laneLength; i ++) {
            everyMWCopy[i] = abs(everyMWCopy[i] - bins[j]);
            
        }
        // print("current bin:");
        // print(bins[j]);
        // print("after subtraction:");
        // Array.print(everyMWCopy);

        Array.getStatistics(everyMWCopy, min, max, mean, stdDev);
        // print("min=" + min);

        for (h = 0; h < laneLength; h ++) {
            if (abs(everyMWCopy[h] - min)<1e-6) {
                binPxValues[j] = h;
                //print(binPxValues[j]);
                //print(j);
            }
        }
    }
    Array.print(binPxValues);
}

calcPxFromBins();


function quantBins() {

}


