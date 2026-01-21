FROM phusion/baseimage:noble-1.0.2

RUN apt-get update && apt-get upgrade -y

# Install Apache
RUN apt-get install -y apache2 git python3-venv gcc libz-dev python3-dev wget sqlite3
COPY 000-default.conf /etc/apache2/sites-available/

# Create runit service for Apache
RUN mkdir -p /etc/service/apache2
COPY apache2.run /etc/service/apache2/run
RUN chmod +x /etc/service/apache2/run && a2enmod cgi

RUN mkdir -p /etc/service/crispor
COPY crispor.run /etc/service/crispor/run
RUN chmod +x /etc/service/crispor/run

# rs3 only works with Python 3.9, so install that now and make venv for it
RUN add-apt-repository ppa:deadsnakes/ppa -y \
 && apt-get install -y python3.9-venv python3.9 python3.9-dev

RUN mkdir -p /data/www \
    && cd /data/www \
    && git clone https://github.com/maximilianh/crisporWebsite.git crispor --depth=1 \
    && cd /data/www/crispor \
    && python3.9 -m venv venv \
    && venv/bin/pip install --upgrade pip setuptools wheel \
    && venv/bin/pip install biopython numpy scikit-learn pandas twobitreader xlwt keras tensorflow h5py rs3 pytabix matplotlib lmdbm PyMySQL \
    && echo 'JOBQUEUEDB = "/data/temp/crisporJobs.db"' > /data/www/crispor/crispor.conf \
    && mkdir -p /data/genomes/ \
    && ln -s /data/genomes /data/www/crispor/genomes \
    && ln -s /data/genomes /data/www/genomes \
    && mkdir -p /data/temp/ \
    && chmod a+rw /data/temp \
    && ln -s /data/temp /data/www/crispor/temp \
    && mkdir /data/www/crispor/log \
    && echo '30 03 * * * root /data/www/crispor/tools/crisprCleanTemp' > /etc/cron.d/crispor

# Clean up APT when done.
# RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Start baseimage-docker's init system.
CMD ["/sbin/my_init"]

