# Gene_Catalog_Rshiny

BAltic Gene Set - BAGS.v1.1

## Deploying to the Serve plattform

### Configuration

The application is configured in the .Renviron.template file. This is initially configured for deployment to the Serve plattform.

### Setup

Upload the dataset to Serve and extract such that the contents reside directly under project-vol. The project-vol directory will be mapped to /data in the container.

Create as custom app type in Serve with the following values:

- Persistent Volume: project-vol
- Path: /data
- Port: 3838
- Image: specify the path to the published image
- User ID: 999


### Viewing logs

    cat /rlogs/app.log

    ls /var/log/shiny-server

## Local deployments

## Build and run a docker image locally

Configure for the local environment in the .Renviron.template file

    R_LOGLEVEL = "DEBUG"
    R_LOGFILE = "../rlogs/app.log"
    DATA_DIR = "../data"

Build the docker image

    docker build -t gene-catalog-rshiny:dev .

Run the container without data folder

    docker run -p 127.0.0.1:3838:3838 gene-catalog-rshiny:dev

Or, run the container with a mounted local data folder

    APP_PATH="PATH-TO-GIT-REPO/Gene_Catalog_Rshiny"

    docker run -p 127.0.0.1:3838:3838 \
    -v $APP_PATH/data:/data \
    --memory=16g --memory-swap=20g \
    gene-catalog-rshiny:dev

Browse to the app at  http://localhost:3838/

## Local development

If planning to run in local environment from R terminal:
Copy and edit the .Renviron file

    cp ./.Renviron.template ./app/.Renviron

If using a local data directory, create the directory and copy a data file there named existing.csv.
This can be used for diagnosing permission problems.

    mkdir data
