Attribute VB_Name = "CynoMathExample"
Option Explicit

' Windows example for Excel/VBA.
' Place cynoMathTools.dll somewhere on the DLL search path,
' or change the declared library path to an absolute file path.

Public Declare PtrSafe Function cyno_BS Lib "cynoMathTools.dll" ( _
    ByVal Tau As Double, _
    ByVal TD As Double, _
    ByVal S As Double, _
    ByVal K As Double, _
    ByVal rD As Double, _
    ByVal rF As Double, _
    ByVal sigma As Double, _
    ByVal callPut As String) As Double

Public Declare PtrSafe Function cyno_Rand1DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long) As Boolean

Public Declare PtrSafe Sub cyno_Rand_Seed Lib "cynoMathTools.dll" ( _
    ByVal seed As Long)

Public Declare PtrSafe Function cyno_Rand1DArrayPlainAntithetic Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long) As Boolean

Public Declare PtrSafe Function cyno_Sobol1DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long) As Boolean

Public Declare PtrSafe Function cyno_SobolNormal1DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long) As Boolean

Public Declare PtrSafe Function cyno_Sobol2DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal dim As Long, _
    ByVal samplesPerFactor As Long) As Boolean

Public Declare PtrSafe Function cyno_SobolNormal2DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal dim As Long, _
    ByVal samplesPerFactor As Long) As Boolean

Public Declare PtrSafe Function cyno_NormalRand1DArrayAntithetic Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_CorrelatedNormalRand1DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef correlationMatrix As Double, _
    ByVal dim As Long, _
    ByRef outArray As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_CorrelatedNormalRand2DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef correlationMatrix As Double, _
    ByVal dim As Long, _
    ByRef outArray As Double, _
    ByVal samplesPerFactor As Long, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_CorrelatedNormalRand1DArrayPlainAntithetic Lib "cynoMathTools.dll" ( _
    ByRef correlationMatrix As Double, _
    ByVal dim As Long, _
    ByRef outArray As Double, _
    ByRef antitheticOutArray As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_CorrelatedNormalRand2DArrayPlainAntithetic Lib "cynoMathTools.dll" ( _
    ByRef correlationMatrix As Double, _
    ByVal dim As Long, _
    ByRef outArray As Double, _
    ByVal samplesPerFactor As Long, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_MatrixMultiplyPlain Lib "cynoMathTools.dll" ( _
    ByRef a As Double, _
    ByVal aRows As Long, _
    ByVal aCols As Long, _
    ByRef b As Double, _
    ByVal bRows As Long, _
    ByVal bCols As Long, _
    ByRef output As Double) As Boolean

Public Declare PtrSafe Function cyno_Sort Lib "cynoMathTools.dll" ( _
    ByRef indexArray As Long, _
    ByRef values As Double, _
    ByVal arraySize As Long) As Long

Public Declare PtrSafe Function cyno_EquityPathPlain Lib "cynoMathTools.dll" ( _
    ByVal spot As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPath As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_EquityPathAntitheticPlain Lib "cynoMathTools.dll" ( _
    ByVal spot As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPath As Double, _
    ByRef antitheticOutPath As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_EquityPathSobolPlain Lib "cynoMathTools.dll" ( _
    ByVal spot As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPath As Double) As Boolean

Public Declare PtrSafe Function cyno_EquityPathsPlain Lib "cynoMathTools.dll" ( _
    ByRef spots As Double, _
    ByVal nStocks As Long, _
    ByRef correlationMatrix As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPaths As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_EquityPathsAntitheticPlain Lib "cynoMathTools.dll" ( _
    ByRef spots As Double, _
    ByVal nStocks As Long, _
    ByRef correlationMatrix As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPaths As Double, _
    ByRef antitheticOutPaths As Double, _
    ByVal randflag As Integer) As Boolean

Public Declare PtrSafe Function cyno_EquityPathsSobolPlain Lib "cynoMathTools.dll" ( _
    ByRef spots As Double, _
    ByVal nStocks As Long, _
    ByRef correlationMatrix As Double, _
    ByRef dt As Double, _
    ByRef drift As Double, _
    ByRef dividend As Double, _
    ByRef sigma As Double, _
    ByVal nSteps As Long, _
    ByRef outPaths As Double) As Boolean

Public Sub RunCynoMathExample()
    Dim callPrice As Double
    callPrice = cyno_BS(1#, 1#, 100#, 95#, 0.03, 0.01, 0.2, "CALL")
    Debug.Print "Black-Scholes call price:", callPrice

    Dim uniforms(0 To 4) As Double
    Dim randomOk As Boolean
    randomOk = cyno_Rand1DArrayPlain(uniforms(0), 5)
    Debug.Print "Random fill ok:", randomOk

    Dim i As Long
    For i = 0 To 4
        Debug.Print "uniform(" & i & ")="; uniforms(i)
    Next i

    Dim sobol(0 To 7) As Double
    Dim sobolOk As Boolean
    sobolOk = cyno_Sobol1DArrayPlain(sobol(0), 8)
    Debug.Print "Sobol fill ok:", sobolOk
    For i = 0 To 7
        Debug.Print "sobol(" & i & ")="; sobol(i)
    Next i

    Dim antitheticNormals(0 To 5) As Double
    Dim antitheticOk As Boolean
    antitheticOk = cyno_NormalRand1DArrayAntithetic(antitheticNormals(0), 6, 1)
    Debug.Print "Antithetic normal fill ok:", antitheticOk
    For i = 0 To 5 Step 2
        Debug.Print "pair(" & (i \ 2) & ")="; antitheticNormals(i); antitheticNormals(i + 1)
    Next i

    Dim left(0 To 3) As Double
    Dim right(0 To 3) As Double
    Dim product(0 To 3) As Double

    left(0) = 1#: left(1) = 2#
    left(2) = 3#: left(3) = 4#

    right(0) = 5#: right(1) = 6#
    right(2) = 7#: right(3) = 8#

    Dim matrixOk As Boolean
    matrixOk = cyno_MatrixMultiplyPlain(left(0), 2, 2, right(0), 2, 2, product(0))
    Debug.Print "Matrix multiply ok:", matrixOk
    Debug.Print "[" & product(0) & " " & product(1) & "]"
    Debug.Print "[" & product(2) & " " & product(3) & "]"
End Sub

Public Sub RunCynoRandomApiExample()
    Dim uniformAntithetic(0 To 5) As Double
    Dim sobolNormal(0 To 5) As Double
    Dim sobol2D(0 To 5) As Double
    Dim sobolNormal2D(0 To 5) As Double
    Dim corr1D(0 To 3) As Double
    Dim correlated(0 To 1) As Double
    Dim correlatedAnti(0 To 1) As Double
    Dim ok As Boolean
    Dim i As Long

    cyno_Rand_Seed 12345

    ok = cyno_Rand1DArrayPlainAntithetic(uniformAntithetic(0), 6)
    Debug.Print "Uniform antithetic fill ok:", ok
    For i = 0 To 5 Step 2
        Debug.Print "u-pair(" & (i \ 2) & ")="; uniformAntithetic(i); uniformAntithetic(i + 1)
    Next i

    ok = cyno_SobolNormal1DArrayPlain(sobolNormal(0), 6)
    Debug.Print "Sobol normal fill ok:", ok
    For i = 0 To 5
        Debug.Print "sobolNormal(" & i & ")="; sobolNormal(i)
    Next i

    ok = cyno_Sobol2DArrayPlain(sobol2D(0), 2, 3)
    Debug.Print "Sobol 2D fill ok:", ok
    Call PrintMatrix("Sobol 2D uniforms", sobol2D, 2, 3)

    ok = cyno_SobolNormal2DArrayPlain(sobolNormal2D(0), 2, 3)
    Debug.Print "Sobol normal 2D fill ok:", ok
    Call PrintMatrix("Sobol 2D normals", sobolNormal2D, 2, 3)

    corr1D(0) = 1#: corr1D(1) = 0.4
    corr1D(2) = 0.4: corr1D(3) = 1#

    ok = cyno_CorrelatedNormalRand1DArrayPlain(corr1D(0), 2, correlated(0), 1)
    Debug.Print "Correlated normal 1D fill ok:", ok
    Debug.Print "correlated(0)="; correlated(0); "correlated(1)="; correlated(1)

    ok = cyno_CorrelatedNormalRand1DArrayPlainAntithetic(corr1D(0), 2, correlated(0), correlatedAnti(0), 1)
    Debug.Print "Correlated normal antithetic fill ok:", ok
    Debug.Print "base vector:"; correlated(0); correlated(1)
    Debug.Print "antithetic vector:"; correlatedAnti(0); correlatedAnti(1)
End Sub

Public Sub RunCynoSortExample()
    Dim values(0 To 3) As Double
    Dim sortedIndex(0 To 3) As Long
    Dim i As Long

    values(0) = 3#
    values(1) = 1#
    values(2) = 2#
    values(3) = -4#

    Call cyno_Sort(sortedIndex(0), values(0), 4)

    Debug.Print "Sorted order from cyno_Sort:"
    For i = 0 To 3
        Debug.Print "position"; i; "index"; sortedIndex(i); "value"; values(sortedIndex(i))
    Next i
End Sub

Public Sub RunCynoNormal2DExample()
    Dim corr(0 To 3) As Double
    Dim normals(0 To 7) As Double
    Dim samplesPerFactor As Long
    Dim factor As Long
    Dim sample As Long
    Dim ok As Boolean

    samplesPerFactor = 4

    ' 2x2 identity matrix gives independent standard normals.
    corr(0) = 1#: corr(1) = 0#
    corr(2) = 0#: corr(3) = 1#

    ok = cyno_CorrelatedNormalRand2DArrayPlain(corr(0), 2, normals(0), samplesPerFactor, 1)
    Debug.Print "2D normal fill ok:", ok

    For factor = 0 To 1
        Debug.Print "factor"; factor
        For sample = 0 To samplesPerFactor - 1
            Debug.Print "  sample"; sample; "value"; normals(factor * samplesPerFactor + sample)
        Next sample
    Next factor
End Sub

Public Sub RunCynoSingleStockPathExample()
    Const nSteps As Long = 4
    Dim dt(0 To nSteps - 1) As Double
    Dim drift(0 To nSteps - 1) As Double
    Dim dividend(0 To nSteps - 1) As Double
    Dim sigma(0 To nSteps - 1) As Double
    Dim path(0 To nSteps) As Double
    Dim stepIndex As Long
    Dim ok As Boolean

    dt(0) = 0.25: dt(1) = 0.25: dt(2) = 0.25: dt(3) = 0.25
    drift(0) = 0.05: drift(1) = 0.05: drift(2) = 0.05: drift(3) = 0.05
    dividend(0) = 0.01: dividend(1) = 0.01: dividend(2) = 0.0125: dividend(3) = 0.0125
    sigma(0) = 0.2: sigma(1) = 0.22: sigma(2) = 0.21: sigma(3) = 0.19

    cyno_Rand_Seed 20260514
    ok = cyno_EquityPathPlain(100#, dt(0), drift(0), dividend(0), sigma(0), nSteps, path(0), 1)

    Debug.Print "Single-stock equity path ok:", ok
    For stepIndex = 0 To nSteps
        Debug.Print "path(" & stepIndex & ")="; path(stepIndex)
    Next stepIndex
End Sub

Public Sub RunCynoThreeStockPathsExample()
    Const nStocks As Long = 3
    Const nSteps As Long = 4
    Dim spots(0 To nStocks - 1) As Double
    Dim corr(0 To nStocks * nStocks - 1) As Double
    Dim dt(0 To nSteps - 1) As Double
    Dim drift(0 To nStocks * nSteps - 1) As Double
    Dim dividend(0 To nStocks * nSteps - 1) As Double
    Dim sigma(0 To nStocks * nSteps - 1) As Double
    Dim outPaths(0 To nStocks * (nSteps + 1) - 1) As Double
    Dim stockIndex As Long
    Dim stepIndex As Long
    Dim baseIndex As Long
    Dim ok As Boolean

    spots(0) = 100#
    spots(1) = 85#
    spots(2) = 120#

    corr(0) = 1#: corr(1) = 0.35: corr(2) = 0.2
    corr(3) = 0.35: corr(4) = 1#: corr(5) = 0.5
    corr(6) = 0.2: corr(7) = 0.5: corr(8) = 1#

    dt(0) = 0.25: dt(1) = 0.25: dt(2) = 0.25: dt(3) = 0.25

    For stockIndex = 0 To nStocks - 1
        For stepIndex = 0 To nSteps - 1
            baseIndex = stockIndex * nSteps + stepIndex
            drift(baseIndex) = 0.04 + 0.005 * stockIndex
            dividend(baseIndex) = 0.01 + 0.002 * stockIndex
            sigma(baseIndex) = 0.18 + 0.03 * stockIndex + 0.005 * stepIndex
        Next stepIndex
    Next stockIndex

    cyno_Rand_Seed 20260514
    ok = cyno_EquityPathsPlain(spots(0), nStocks, corr(0), dt(0), drift(0), dividend(0), sigma(0), nSteps, outPaths(0), 1)

    Debug.Print "Three-stock equity paths ok:", ok
    For stockIndex = 0 To nStocks - 1
        Debug.Print "stock"; stockIndex
        For stepIndex = 0 To nSteps
            Debug.Print "  path(" & stepIndex & ")="; outPaths(stockIndex * (nSteps + 1) + stepIndex)
        Next stepIndex
    Next stockIndex
End Sub

Private Sub PrintMatrix(ByVal title As String, ByRef matrix() As Double, ByVal rows As Long, ByVal cols As Long)
    Dim rowIndex As Long
    Dim colIndex As Long
    Dim lineText As String

    Debug.Print title
    For rowIndex = 0 To rows - 1
        lineText = ""
        For colIndex = 0 To cols - 1
            lineText = lineText & Format$(matrix(rowIndex * cols + colIndex), "0.000000") & " "
        Next colIndex
        Debug.Print Trim$(lineText)
    Next rowIndex
End Sub
