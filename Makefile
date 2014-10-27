TARGET = gsmf
SRCDIR = src
all p prof d debug:
	@cd $(SRCDIR) && make $(MAKECMDGOALS)
	@-[ -e $(SRCDIR)/$(TARGET) ] && cp $(SRCDIR)/$(TARGET) . ||:

c clean:
	@cd $(SRCDIR) && make $(MAKECMDGOALS)
	@-[ -e $(SRCDIR)/$(TARGET) ] && cp $(SRCDIR)/$(TARGET) . ||:
	@-[ -e $(TARGET) ] && rm -v $(TARGET) ||:

