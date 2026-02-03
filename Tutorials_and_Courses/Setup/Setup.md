# Setup
## Summary
***

This tutorial will teach you how to connect to a server, transfer files between your local computer and the server, basic commands, and running an interactive session. 

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
  <img src=https://github.com/felipehcoutinho/ICM_Code/blob/main/Tutorials_and_Courses/Figures/putty_config.png" width="400" height="400" alt="Putty config example"/>
</p>

- Run WinSCP. Click New site, and fill the username, address and password fields. Click save to store this for future use.

<p align="center">
  <img src="https://github.com/felipehcoutinho/ICM_Code/blob/main/Tutorials_and_Courses/Figures/winscp_config.png" width="400" height="600" alt="WiNSCp config example"/>
</p>

Linux and Mac users:
- Open the terminal an run the command below replacing 'username' and 'address' accordingly to the server you are connecting to:
`ssh username@address`

## Connecting to MARBITS from outside the ICM Network:

https://marbits.icm.csic.es/documentation/


## Introduction to shell commands

https://www.geeksforgeeks.org/linux-unix/basic-shell-commands-in-linux/

## Bash tutorial

Recommended studying the Introduction and Scripting modules

https://www.w3schools.com/bash/index.php 

## Introduction to Slurm

https://blog.ronin.cloud/slurm-intro/

### Starting a slurm interactive session

Start a session requesting 4 CPUs, 5 gigabytes of RAM, and running for a maximum time of 3 hours, 33 minutes, and 03 seconds:

`srun --threads 4 --mem 5G --time 03:33:03 --pty /bin/bash`

End the session with the following command

`exit`



