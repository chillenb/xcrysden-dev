default:

distclean:
	for dir in bin etc include lib lib64 share usr man; do \
		if test -d $$dir ; then rm -rf $$dir ; fi; \
	done
	cd src; $(MAKE) clean
