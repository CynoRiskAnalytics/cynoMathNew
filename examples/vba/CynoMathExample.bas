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

Public Declare PtrSafe Function cyno_Sobol1DArrayPlain Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long) As Boolean

Public Declare PtrSafe Function cyno_NormalRand1DArrayAntithetic Lib "cynoMathTools.dll" ( _
    ByRef outArray As Double, _
    ByVal size As Long, _
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
