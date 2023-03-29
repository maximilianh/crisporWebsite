How to clean a Crispor server

    find /data/www/temp -type f -mtime +30 | egrep '(effScores.tab)|(primer3)|(bed.gz)|(.db$)' | sudo xargs rm
