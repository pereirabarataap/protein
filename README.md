## Usage example

<code>
  This function performs the following, in sequence:
        1. Compute correlation matrix of protein-protein from protein_df
        2. Hierarchical clustering of proteins using each correlation-row as feature vector 
            2.1 Ward-linkage is used to compute distances between proteins
            2.2 Optimal number of clusters is selected on lowest Davies Bouldin score
        3. Venn diagram of each cluster showing GO terms of respective cluster proteins
        4. Counts number of proteins in each cluster which contain cluster-specific GO terms
        
    protein_df ---> pandas.DataFrame
                    Column names must be protein accessions: str like "P05067"
                    Row values represent measures of respective protein: float
                    
        |   |     P98160 |   O00468-6 |   Q15149-4 |     P08238 |     P06733 |
        |:--|-----------:|-----------:|-----------:|-----------:|-----------:|
        | 0 | 1.5315e+09 | 2.5305e+09 | 9.382e+08  | 2.9324e+09 | 4.9473e+09 |
        | 1 | 1.8858e+09 | 3.1737e+09 | 8.2028e+08 | 3.4993e+09 | 3.6529e+09 |
        | 2 | 1.1927e+09 | 2.5913e+09 | 6.1717e+08 | 2.7219e+09 | 7.6858e+09 |
        | 3 | 1.4522e+09 | 2.6129e+09 | 9.3064e+08 | 3.2274e+09 | 3.6938e+09 |
        | 4 | 1.6498e+09 | 2.3794e+09 | 6.6648e+08 | 1.9489e+09 | 4.2946e+09 |
</code>

![image](https://user-images.githubusercontent.com/15198092/164678427-f07ceaff-27e1-417d-bfe4-d4e19e821fb3.png)

