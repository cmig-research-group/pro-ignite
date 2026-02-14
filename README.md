## How to deploy Pro-IGNITE software
### Prerequisites
Linux: These instructions assume your system is running Ubuntu. Other Linux distrubutions will likely work with minor changes (e.g., apt -> dnf for RHEL)

### Clone this repo
```
git clone https://github.com/cmig-research-group/pro-ignite.git
```

### Install Docker
```
sudo apt install docker.io
```

### Pull Docker images
The core Pro-IGNITE software is packaged as two Docker images, hosted on GitHub [here](https://github.com/orgs/cmig-research-group/packages). To pull them onto your local machine, you need a GitHub account with an access token that allows for package access; see documentation [here](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry).
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
PATH_WRITE_TARGET = '/home/orthanc/pro_ignite'
PATH_PROIGNITE_TMP = '/home/tmp_pro_ignite'
PATH_CONFIG_DIR = '/home/ccconlin/work/pro-ignite'
PATH_LOGS = '/home/orthanc/pro_ignite_logs'
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
