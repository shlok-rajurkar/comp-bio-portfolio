function main(workbook: ExcelScript.Workbook) {
  
  const worksheet = workbook.getWorksheet("Sheet1");
  if (!worksheet) {
    console.log("Worksheet 'Sheet1' not found")
    return
  }

  const dataRange = worksheet.getUsedRange(true);
  if (!dataRange) {
    console.log("No data found")
    return
  }

  const values = dataRange.getValues();

  const outputData: (string | number)[][] = [];

  outputData.push([
    "Sample ID",
    "Center",
    "Area"
  ])

  const inputColA = 0;
  const inputColB = 1;
  const inputColD = 3;

  const lowerBound = 14.5;
  const upperBound = 18;

  let sampleID = 0;
  let numMZCurves = 0;
  let mZCurveArea: number | string = 0;
  let mZCurveCenter: number | string = 0;
  for (let i = 0; i < values.length; i++) {

    const row = values[i];

    if (row.length <= inputColD) {
      continue;
    }

    if (row[1] == "Center") {
        sampleID++;
        continue;
    }

    if (row[1] > lowerBound && row[1] < upperBound) {
        mZCurveCenter = row[1];
        mZCurveArea = row[3];
        numMZCurves++;
    }

    if (row[0] == "") {
        if (sampleID > 0) {
            if (numMZCurves === 1) {
                outputData.push([
                sampleID,
                mZCurveCenter, 
                mZCurveArea
                ])
            } else if (numMZCurves === 0) {
                outputData.push([
                sampleID,
                "none",
                "none"
                ])
            } else {
                outputData.push([
                sampleID,
                "excess",
                "excess"
                ])
            }
        }
        numMZCurves = 0;
        mZCurveCenter = "none";
        mZCurveArea = "none"
    }



  }

  if (outputData.length > 0) {
    worksheet.getRange("J:L").clear();

    const outputRange = worksheet.getRangeByIndexes(
      0, 9, outputData.length, 3);
      outputRange.setValues(outputData);

  }


}