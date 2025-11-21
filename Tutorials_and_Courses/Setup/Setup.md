# Setup
## Summary
***

This tutorial will teach you how to connect to a server to run commands and transfer files between your local computer and the server. 

***

## Software
Download and install the following program:

- [VSCode](https://code.visualstudio.com/download)

Note: We will use VSCode to make it easier to write and edit code. It is VERY important to save all the commands you run (with comments explaining each step) so that later you can check exactly how the data was treated and how the results were generated. 

Windows users only:
Also download and install the following programs that you will use to connect and transfer files to the server:
- [WinSCP](https://winscp.net/eng/download.php)
- [Putty](https://putty.org/index.html)

## Connecting to MARBITS from inside the ICM Network:
Windows users:
- Run putty, fill the Host name (or IP address) field with: username@address  replacing 'username' and 'address' accordingly to the server you are connecting to. Fill the "saved sessions" field with an easy to remember name for the server. Click save to store this for future use.

<p align="center">
  <img src="https://github.com/felipehcoutinho/virathon/blob/main/Virathon_Logo.png" width="400" height="400" alt="Virathon logo generated with DALLE"/>
</p>

- Run WinSCP. Click New site, and fill the username, address and password fields. Click save to store this for future use.

<p align="center">
  <img src="https://github.com/felipehcoutinho/virathon/blob/main/Virathon_Logo.png" width="400" height="400" alt="Virathon logo generated with DALLE"/>
</p>

Linux and Mac users:
- Open the terminal an run the command below replacing 'username' and 'address' accordingly to the server you are connecting to:
`ssh username@address`

## Connecting to MARBITS from outside the ICM Network:

https://marbits.icm.csic.es/documentation/

