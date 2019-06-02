# compara-deep-learning
Using Deep Learning techniques to enhance orthology calls

This project aims to apply machine learning algorithms like Deep Learning Neural Networks to validate the homologies predicted with our method in addition to infer new ones based on other properties of the data that are currently not being considered (such as local synteny, divergence rates, etc).

Run the `ftpg.py` directly to get the links of the files in the ftp server to be downloaded. Enter `y` when prompted to download the files to the files. The directories to which the files are downloaded can be modified in the script.

The command line arguments are as follows:<br/>
`-f` read the links from the file and download it<br/>
`-d` read the files in the given directory<br/>
`-nd` ignore this argument(do nothing)<br/>
`-l` download file from the link. Works <b>only for homology</b> file<br/>
`-r` or `-d`  to specify to download and read or just download data<br/>
`-test` or `-run` to test the files and the code or run it.<br/>

<b>For Example:</b><br/>
  `python main.py -f link.txt -d data_homology -r -test` will download the files given by the links in the file `link.txt` to the  `data` directory and read the downloaded files.It will also read all the files in the `data_homology` directory. The last argument shows that the  data will be downloaded as well as read. All the files will be read and selective records will be processed.
  
The sample output for `-test` should match this at the end:<br/>
Synteny Matrices are created successfully<br/>
 3<br/>
 3<br/>
(3, 21)<br/>
(3, 21)<br/>
(3,)<br/>
(3,)<br/>
(3, 7, 7, 2)<br/>
(3, 7, 7, 2)<br/>
(3,)<br/>
(3,)<br/>
(3,)<br/>
Data Saved Successfully to processed :)<br/>

At this moment the model has been designed to work with <b>only one</b> homology database. Functionality will be updated during the course of the project.


