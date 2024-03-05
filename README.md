# Gene_Catalog_Rshiny
Fork of lfdelzam/Gene_Catalog_Rshiny

For testing variations of building and deploying this Shiny app.

## Development

If planning to run in local environment from R terminal:
Copy and edit the .Renviron file

    cp ./.Renviron.template ./app/.Renviron

If using a local data directory, create the directory and copy a data file there named existing.csv:

    mkdir data

## Build a docker image locally

    docker build -t gene-catalog-rshiny:dev .

Run the container

    docker run -p 127.0.0.1:3838:3838 gene-catalog-rshiny:dev

Browse to the app at  http://localhost:3838/

## View the log files

    cat /rlogs/app.log

    ls /var/log/shiny-server
