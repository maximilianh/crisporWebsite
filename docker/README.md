The CRISPOR docker container is based on Ubuntu 26 (Phusion, for docker) and
includes Apache, Sqlite3, Python 3.9, a ton of Python dependencies for the
efficiency scoring algorithms (keras, numpy, numba, scikit-learn, etc), a
cronjob to clean up the temp files and CRISPOR itself. The software is
installed under /data/www/crispor, which is also the htdocs directory for
Apache. The scoring daemon is started on boot under /etc/system.d/crispor by
initd. 

To download and start the container and map the port 8080 on your machine to the container:

   docker run -p 8080:80 --name crispor-container maximilianh/crispor /sbin/my_ini

You should then be able to access the container via http://localhost from the machine. There is no genome yet.

To download a genome from crispor.gi.ucsc.edu into the container:

   docker exec -it crispor-container /data/www/crispor/tools/crisporDownloadGenome hg38

To add a new genome to the container:

   docker exec -it crispor-container /data/www/crispor/tools/crisporDownloadGenome hg38

To inspect the docker container:

   docker exec -it crispor-container333 /bin/bash

I run this command to build the container and push it as a multi-architecture build:

   docker buildx build . --platform linux/amd64,linux/arm64 -t maximilianh/crispor:latest --push

Then I tag the container with the version, since it's a multi-arch container, I cannot use the tag comment, but need to use buildx:

   docker buildx imagetools create --tag maximilianh/crispor:v5.2 maximilianh/crispor:latest
