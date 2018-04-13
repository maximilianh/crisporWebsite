git:
	sed -i '1 s|^.*$$|#!/usr/bin/env python2.7|' crispor.cgi
dev:
	sed -i '1 s|^.*$$|#!/cluster/software/bin/python2.7|' crispor.cgi

commit:
	make git; git add crispor.cgi; git commit; make dev;
