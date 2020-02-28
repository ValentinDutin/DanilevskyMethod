import java.util.ArrayList;
import java.util.List;

public class DanilevskyMethod {
    private double symmetrickMatr[][];
    private double matrA [][];
    private double matrM [][];
    private double inverseM [][];
    private double vectorY[];
    private double matrS[][];
    private int n;
    private List<double[]> eigenVectors;
    private double lambda [] = {0.1910090105395, 0.383558316493, 0.597703453394, 0.879661472553, 1.144756197021};
    private List<double[]> discrepancy;
    private double[] polinom;


    public DanilevskyMethod(double matrA[][]){
        n = matrA.length;
        symmetrickMatr = new double[n][n];
        matrS = new double[n][n];
        this.matrA = new double[n][n];
        matrM = new double[n][n];
        inverseM = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                this.matrA[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    this.matrA[i][j] += matrA[k][i] * matrA[k][j];
                }
                symmetrickMatr[i][j] = this.matrA[i][j];
            }
        }
    }

    private void createMatrM(int count){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == n - 2 - count) {
                    if (j == i) {
                        matrM[i][j] = 1 / matrA[i+1][i];
                    } else {
                        matrM[i][j] = -matrA[i+1][j] / matrA[i+1][i];
                    }
                }
                else if(i == j){
                    matrM[i][i] = 1;
                }
                else{
                    matrM[i][j] = 0;
                }
                if(count == 0){
                    matrS[i][j] = matrM[i][j];
                }
            }
        }
    }
    private void createInverseM(int count){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == n - 2 - count) {
                    inverseM[i][j] = matrA[i+1][j];
                }
                else if(i == j){
                    inverseM[i][i] = 1;
                }
                else{
                    inverseM[i][j] = 0;
                }
            }
        }
    }
    private double[][] multiply(double matrA[][], double matrB[][]){
        double result[][] = new double[n][n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = 0;
                for (int k = 0; k < n; k++)
                    result[i][j] += matrA[i][k] * matrB[k][j];
            }
        }
        return result;
    }
    private double[] multiply(double[][] matr, double[] vector){
        double result[] = new double[n];
        for(int i=0; i<n; i++)
        {
            result[i]=0;
            for(int j=0; j<n; j++)
            {
                result[i] += matr[i][j]*vector[j];
            }
        }
        return result;
    }
    private double[] multiply(double[] vector, double lambda){
        for(int i = 0; i < n; i++){
            vector[i] *= lambda;
        }
        return vector;
    }
    private double[] minus(double[] vectorA, double[] vectorB){
        for(int i = 0; i < n; i++){
            vectorA[i] -= vectorB[i];
        }
        return vectorA;
    }

    public void danilevskyMethod(){
        for(int step = 0; step < n-1; step++){
            createMatrM(step);
            createInverseM(step);
            matrA = multiply(multiply(inverseM, matrA), matrM);
            if(step > 0){
                matrS = multiply(matrS, matrM);
            }
        }
        polinom = new double[n + 1];
        polinom[0] = 1;
        for(int i = 1; i < n+1; i++){
            polinom[i] = -matrA[0][i-1];
        }
    }

    public void printA(){
        for (double[] row: matrA){
            for(double item: row){
                System.out.format("%25s", item + "    ");
            }
            System.out.println();
        }
        System.out.println();
    }
    public void createEigenVectors(){
        eigenVectors = new ArrayList<>();
        for(int i = 0; i < n; i++){
            eigenVectors.add(createEigenVector(lambda[i]));
        }
    }

    public void printDiscrepancy(){
        double discrepancy[] = new double[n];
        for(int i = 0; i < n; i++){
            discrepancy[i] = 0;
            for(int j = 0; j < n; j++){
                discrepancy[i] += Math.pow(lambda[i], n-j) * polinom[j];
            }
            discrepancy[i] += polinom[n];
            System.out.println("discrepancy for lambda = " + lambda[i]);
            System.out.format("%25s", discrepancy[i] + "\n");
        }
    }

    private double[] createEigenVector(double lambda){
        vectorY = new double[n];
        for(int i = 0; i < n; i++){
            vectorY[i] = Math.pow(lambda, n-1-i);
        }
        return multiply(matrS, vectorY);
    }
    public void printEigenVectors(){
        int count = 0;
        for(double[] vector: eigenVectors){
            System.out.println("Eigen vector for lambda = " + lambda[count]);
            for(double item : vector){
                System.out.format("%25s", item + "    ");
            }
            System.out.println();
            count++;
        }
    }

    private double[] createVectorsDiscrepancy(double[] eigenVector, double lambda){
        return minus(multiply(symmetrickMatr, eigenVector), multiply(eigenVector, lambda));
    }
    public void createVectorsDiscrepancy(){
        discrepancy = new ArrayList<>();
        for(int i = 0; i < n; i++){
            discrepancy.add(createVectorsDiscrepancy(eigenVectors.get(i), lambda[i]));
        }
    }
    public void printVectorsDiscrepancy(){
        int count = 0;
        for(double[] vector: discrepancy){
            System.out.println("discrepancy for lambda = " + lambda[count]);
            for(double item : vector){
                System.out.format("%25s", item + "    ");
            }
            System.out.println();
            count++;
        }
    }
}
