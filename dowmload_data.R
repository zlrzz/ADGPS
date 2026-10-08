# sample information download form the website below
# https://www.ncbi.nlm.nih.gov/geo/browse/
# choose Samples and search the Series number of GEO, then Export file to local destination
# GSE5281 -----------------------------------------------------------------

if(dir.exists("data/GSE5281") == F) {
  dir.create("data/GSE5281")
}

if(file.exists("data/GSE5281/GSE5281_series_matrix.txt.gz") == F) {
  downloader::download(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5281/matrix/GSE5281_series_matrix.txt.gz",
                       destfile = "data/GSE5281/GSE5281_series_matrix.txt.gz")
}

if(file.exists("data/GSE5281/sample.csv") == F) {
  print("Please dowmload sample.csv for this GEO session.")
}

# GSE48350 ----------------------------------------------------------------

if(dir.exists("data/GSE48350") == F) {
  dir.create("data/GSE48350")
}

if(file.exists("data/GSE48350/GSE48350_series_matrix.txt.gz") == F) {
  downloader::download(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/GSE48350_series_matrix.txt.gz",
                       destfile = "data/GSE48350/GSE48350_series_matrix.txt.gz")
}

if(file.exists("data/GSE48350/sample.csv") == F) {
  print("Please dowmload sample.csv for this GEO session.")
}

# GSE159699 ---------------------------------------------------------------

if(dir.exists("data/GSE159699") == F) {
  dir.create("data/GSE159699")
}

if(file.exists("data/GSE159699/GSE159699_summary_count.star.txt.gz") == F) {
  downloader::download(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159699&format=file&file=GSE159699%5Fsummary%5Fcount%2Estar%2Etxt%2Egz",
                       destfile = "data/GSE159699/GSE159699_summary_count.star.txt.gz")}

if(file.exists("data/GSE159699/sample.csv") == F) {
  print("Please dowmload sample.csv for this GEO session.")
}

# GSE104704 ---------------------------------------------------------------

if(dir.exists("data/GSE104704") == F) {
  dir.create("data/GSE104704")
}

if(file.exists("data/GSE104704/GSE104704_RNA-Seq_Table.txt.gz") == F) {
  downloader::download(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE104704&format=file&file=GSE104704%5FRNA%2DSeq%5FTable%2Etxt%2Egz",
                       destfile = "data/GSE104704/GSE104704_RNA-Seq_Table.txt.gz")
}

if(file.exists("data/GSE104704/sample.csv") == F) {
  print("Please dowmload sample.csv for this GEO session.")
}

# GSE33000 ----------------------------------------------------------------

if(dir.exists("data/GSE33000") == F) {
  dir.create("data/GSE33000")
}

if(file.exists("data/GSE33000/GSE33000_raw_data.txt.gz") == F) {
  downloader::download(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE33000&format=file&file=GSE33000%5Fraw%5Fdata%2Etxt%2Egz",
                       destfile = "data/GSE33000/GSE33000_raw_data.txt.gz")
}

if(file.exists("data/GSE33000/GSE33000_series_matrix.txt.gz") == F) {
  downloader::download(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33000/matrix/GSE33000_series_matrix.txt.gz",
                       destfile = "data/GSE33000/GSE33000_series_matrix.txt.gz")
}

if(file.exists("data/GSE33000/sample.csv") == F) {
  print("Please dowmload sample.csv for this GEO session.")
}

# end ---------------------------------------------------------------------

print("Had downloaded all the files you write in this script.")
