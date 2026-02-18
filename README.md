## Overview of Pro-IGNITE software architecture
![](/software_architecture.png)

## How to deploy Pro-IGNITE software
### Prerequisites
Linux: These instructions assume your system is running Ubuntu. Other Linux distrubutions will likely work with minor changes (e.g., apt -> dnf for RHEL)

### Clone this repo
```
git clone https://github.com/cmig-research-group/pro-ignite.git
```

### Configure Pro-IGNITE processing
There are 3 configuration files for different parts of the software:

**`config_app.m`** -> Settings for the gestalt app

**`config_routing.m`** -> Defines DICOM destination(s) where output images should be sent

`config_rsi.m` -> RSI processing options (in general, this file should be left alone)

### Install Docker
```
sudo apt install docker.io
```

### Pull Docker images
The core Pro-IGNITE software is packaged as two Docker images, hosted on GitHub [here](https://github.com/orgs/cmig-research-group/packages). Install them via `docker pull`:
```
docker pull ghcr.io/cmig-research-group/autoseg_prostate:latest
docker pull ghcr.io/cmig-research-group/pro_ignite:latest
```

### Install Orthanc DICOM server
```
sudo apt install orthanc
```

### Configure Orthanc
1. Orthanc will be managed as a systemd service. Stop Orthanc while it's being configured:
```
sudo service orthanc stop
```

2. Install the systemd unit file for Orthanc from this repo: [pro-ignite/orthanc/orthanc.service](https://github.com/cmig-research-group/pro-ignite/blob/main/orthanc/orthanc.service)


3. Install the Orthanc configuration file from this repo by copying it into /etc/orthanc/
```
sudo cp pro-ignite/orthanc/orthanc.json /etc/orthanc/
```

4. Make some directories that will be used by Orthanc (assuming default configuration provided by this repo)
```
sudo mkdir -p /home/orthanc/tmp
sudo mkdir -p /home/orthanc/scripts
```

5. Configure paths at the top of [pro_ignite.lua](https://github.com/cmig-research-group/pro-ignite/blob/main/orthanc/lua/pro_ignite.lua)
```
-- Change this path! It needs to point to the directory where config_app.m is stored.
PATH_CONFIG_DIR = '/home/ccconlin/work/pro-ignite'

-- Only change these paths if you aren't happy with the default configuration
PATH_WRITE_TARGET = '/home/orthanc/pro_ignite'
PATH_LOGS = '/home/orthanc/pro_ignite_logs'
PATH_PROIGNITE_TMP = '/home/tmp_pro_ignite'
```

6. Install [pro_ignite.lua](https://github.com/cmig-research-group/pro-ignite/blob/main/orthanc/lua/pro_ignite.lua) in the directory for Orthanc scripts (assuming default configuration provided by this repo)
```
sudo cp pro-ignite/orthanc/lua/pro_ignite.lua /home/orthanc/scripts/
```

7. Restart Orthanc
```
sudo systemctl daemon-reload
sudo service orthanc restart
```
