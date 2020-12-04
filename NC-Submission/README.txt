###Instructions for using custom code for Smith et al., Nature Communications submission NCOMMS-20-45565
###Note: text following $ in these instructions is meant to be entered into the terminal/command prompt exactly (case sensitive)
###This code relies on outputs from the ClearMap package. The ClearMap documentation with installation instructions is included here,
###and code can be found at www.github.com/christophkirst/clearmap

###System requirements:
Tested on Linux Ubuntu 18.04LTS, but should work on Mac OSX & Windows.
Python libraries specified in requirements.txt

###Installation instructions:
1) Open a terminal/command prompt inside the NatComm_Supplement folder.
2) sudo apt install python3-pip
2) Install dependencies (approximately 2-3 minutes depending on internet connection speed):
    $ python3 pip install -r requirements.txt
3) sudo chmod +x clearMapSubregionParser.py

###Run the code (expected time 5-15 seconds):
    $ python3 clearMapSubregionParser.py --directory ./  --hemisphere 'left' --samples IA1_LB --save True

###Expected Text Output:

Validated correct split
    mouse  aDLS  mDLS  pDLS  aDMS  mDMS  pDMS
0  IA1_LB   151   283   157   149   410   172

###Expected CSV file output: 
./IA1_LB_Striatum_Subregion_Counts_left.csv
