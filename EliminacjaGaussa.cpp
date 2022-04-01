#include #include #include
using namespace std;

// epsilon do przybliżania wartości do 0
const double epsilon = 0.000000001;

// funkcja wypełniająca macierz danymi
void fillMatrix(double **&matrix, double data[], int n, int m) {
int globalI = 0;
for (int i = 0; i < n; i++) {
for (int j = 0; j < m; j++) {
matrix[i][j] = data[globalI++];
}
}
}

// funkcja przydzielająca pamięć macierzy
void allocateMatrix(double **&matrix, int n) {
for (int i = 0; i < n; i++) {
matrix[i] = new double[n + 1]; // [n + 1] - dodatkowa kolumna dla macierzy B
}
}

// funkcja usuwająca macierz
void deleteMatrix(double **&matrix, int n) {
for (int i = 0; i < n; i++) {
delete [] matrix[i];
}
delete [] matrix;
matrix = nullptr;
}

// funkcja wyświetlająca macierz
void showMatrix(double **matrix, int n, int m) {
for (int i = 0; i < n; i++) {
for (int j = 0; j < m; j++) {
cout << setw(10) << matrix[i][j];
}
cout << endl;
}
}

// funckja sprawdza czy liczba = 0, jeśli tak, to zamienia tę liczbę na 0
bool checkIfZeroAndReplace(double &a) {
if (abs(a) < epsilon) {
a = 0;
return true;
}
return false;
}

// liczenie metodą gaussa - wariant 1 (bez zamiany kolumn i wierszy)
// bool, ponieważ na głównej przekątnej może pojawić się 0
// gdy tak się stanie, funkcja zwróci fałsz
bool gauss(double **&matrix, double *X, int n) {
for (int k = 0; k < n-1; k++) { // n-1 kroków k
for (int i = k+1; i < n; i++) { // i = k+1, ..., n-1 (numerowanie od 0)
if (checkIfZeroAndReplace(matrix[k][k]))
return false; // jeœli matrixAB[k][k] = 0, zakończ algorytm
double p = matrix[i][k] / matrix[k][k];
for (int j = 0; j < n+1; j++) { // przejście przez cały wiersz
matrix[i][j] -= p * matrix[k][j];
}
}
}

// postępowanie odwrotne gaussa
X[n-1] = matrix[n-1][n] / matrix[n-1][n-1];
for (int i = n-2; i >= 0; i--) {
double sum = 0;
for (int s = i+1; s < n; s++)
sum += matrix[i][s] * X[s];
X[i] = (matrix[i][n] - sum) / matrix[i][i];
}
if (abs(matrix[n-1][n-1]) <= epsilon)
return false;
return true;
}

// liczenie metodą gaussa - wariant 2 (zamiana wierszy)
bool gauss2(double **&matrix, double *X, int n) {
for (int k = 0; k < n-1; k++) { // n-1 kroków k
int maxTempIndex = k;
for (int i = k; i < n; i++) { // wyszukiwanie największej wartości w kolumnie
if (abs(matrix[i][k]) > abs(matrix[maxTempIndex][k])) {
maxTempIndex = i;
}
}
if(maxTempIndex != k) { // zamień wiersze jeśli indeksy się różnią
swap(matrix[k], matrix[maxTempIndex]);
}
for (int i = k+1; i < n; i++) { // i = k+1, ..., n-1 (numerowanie od 0)
if (checkIfZeroAndReplace(matrix[k][k]))
return false; // jeœli matrixAB[k][k] = 0, zakończ algorytm
double p = matrix[i][k] / matrix[k][k];
for (int j = 0; j < n+1; j++) { // przejście przez cały wiersz
matrix[i][j] -= p * matrix[k][j];
}
}
}

// postępowanie odwrotne gaussa
X[n-1] = matrix[n-1][n] / matrix[n-1][n-1];
for (int i = n-2; i >= 0; i--) {
double sum = 0;
for (int s = i+1; s < n; s++)
sum += matrix[i][s] * X[s];
X[i] = (matrix[i][n] - sum) / matrix[i][i];
}
if (abs(matrix[n-1][n-1]) <= epsilon)
return false;
return true;
}

// funkcja zamieniająca kolumny
void swapColumns(double**& matrix, int n, int col1, int col2) {
for (int i = 0; i < n; i++) {
double temp = matrix[i][col1];
matrix[i][col1] = matrix[i][col2];
matrix[i][col2] = temp;
}
}

// liczenie metodą gaussa - wariant 3 (pełny)
bool gauss3(double **&matrix, double *X, int n) {
int *Xindex = new int[n]; // tablica przechowująca indeksy dla X
for (int i = 0; i < n; i++) { // wypełnienie indeksami tablicy indeksów
Xindex[i] = i;
}

for (int k = 0; k < n-1; k++) { // n-1 kroków k
int maxTempIndexI = k;
int maxTempIndexJ = k;
for (int i = k; i < n; i++) { // wyszukiwanie największej wartości w macierzy pomniejszonej
for (int j = k; j < n; j++) {
if (abs(matrix[i][j]) > abs(matrix[maxTempIndexI][maxTempIndexJ])) {
maxTempIndexI = i;
maxTempIndexJ = j;
}
}
}
if(maxTempIndexI != k) { // zamień wiersze jeśli indeksy się różnią
swap(matrix[k], matrix[maxTempIndexI]);
}
if (maxTempIndexJ != k) { // zamień kolumny jeśli indeksy się różnią
swapColumns(matrix, n, k, maxTempIndexJ);
swap(Xindex[k], Xindex[maxTempIndexJ]); // zamiana indeksow wynikowych
}

for (int i = k+1; i < n; i++) { // i = k+1, ..., n-1 (numerowanie od 0)
if (checkIfZeroAndReplace(matrix[k][k]))
return false; // jeœli matrixAB[k][k] = 0, zakończ algorytm
double p = matrix[i][k] / matrix[k][k];
for (int j = 0; j < n+1; j++) { // przejście przez cały wiersz
matrix[i][j] -= p * matrix[k][j];
}
}
}

// postępowanie odwrotne gaussa
X[n-1] = matrix[n-1][n] / matrix[n-1][n-1];
for (int i = n-2; i >= 0; i--) {
double sum = 0;
for (int s = i+1; s < n; s++)
sum += matrix[i][s] * X[s];
X[i] = (matrix[i][n] - sum) / matrix[i][i];
}

// ustawienie odpowiednich indeksów X
double *XCopy = new double[n];
for (int i = 0; i < n; i++) { // kopiowanie tablicy X
XCopy[i] = X[i];
}
for (int i = 0; i < n; i++) {
X[Xindex[i]] = XCopy[i];
}
if (abs(matrix[n-1][n-1]) <= epsilon)
return false;
return true;
}

// funkcja wyświetla wynik
void showResult(double **matrix, double *X, int n) {
cout << endl;
cout << "Utworzona macierz trojkatna:\n";
showMatrix(matrix, n, n+1);
cout << endl;
cout << "X[1, ..., " << n << "]:\t";
for (int j = 0; j < n; j++) {
checkIfZeroAndReplace(X[j]);
cout << X[j] << "\t";
}
}

int main() {
int choice = 0;
while (choice != 1 && choice != 2) {
cout << "1. Dane \"na sztywno\"" << endl;
cout << "2. Dane z klawiatury\n>> ";
cin >> choice;
}

int method = 0;
while (method != 1 && method != 2 && method != 3) {
cout << "\n\n1. Gauss podstawowy" << endl;
cout << "2. Gauss z zamiana wierszy\n";
cout << "3. Gauss pelny\n>> ";
cin >> method;
}

int n = 0; // wymiar macierzy A
double** matrixAB; // macierz [A|B]
double* X; // macierz wynikowa X

if (choice == 1) {
n = 4;
matrixAB = new double*[n];
allocateMatrix(matrixAB, n);
X = new double[n];
double data[4*5] = {
1, 2, -1, 2, 0,
1, 0, -2, 4, 4,
0, -3, 1.5, 7, 0,
0, -1, 1, 6, -1
};
fillMatrix(matrixAB, data, n, n+1);
} else {
while (n <= 0) {
cout << "Podaj wielkosc macierzy A - n ( > 0): ";
cin >> n;
}
matrixAB = new double*[n];
allocateMatrix(matrixAB, n);
X = new double[n];
cout << "Wprowadzanie macierzy [A|B]\n";
cout << "Macierz A:\n";
for (int i = 0; i < n; i++) {
for (int j = 0; j < n; j++) {
cout << "A[" << i << "][" << j << "]: ";
cin >> matrixAB[i][j];
}
}
cout << "Macierz B:\n";
for (int i = 0; i < n; i++) {
cout << "B[" << i << "]: ";
cin >> matrixAB[i][n];
}
}

cout << "Wprowadzona macierz: \n";
showMatrix(matrixAB, n, n+1);

switch (method) {
case 1:
if (gauss(matrixAB, X, n)) {
showResult(matrixAB, X, n);
} else {
showMatrix(matrixAB, n, n+1);
cout << "Element na glownej przekatnej wynosi 0\n";
}
break;
case 2:
if (gauss2(matrixAB, X, n)) {
showResult(matrixAB, X, n);
} else {
showMatrix(matrixAB, n, n+1);
cout << "Element na glownej przekatnej wynosi 0\n";
}
break;
case 3:
if (gauss3(matrixAB, X, n)) {
showResult(matrixAB, X, n);
} else {
showMatrix(matrixAB, n, n+1);
cout << "Element na glownej przekatnej wynosi 0\n";
}
break;
}

deleteMatrix(matrixAB, n);
delete [] X;

return 0;
}