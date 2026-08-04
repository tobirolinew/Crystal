FFHASH_REPO = https://raw.githubusercontent.com/jannisteunissen/ffhash
FFHASH_BRANCH = master
FFHASH_FILE = ffhash_inc.f90
LICENSE_FILE = LICENCE.txt

.PHONY: all clean

all:
	@if [ -f "$(FFHASH_FILE)" ]; then \
		msg1="!!!! Warning: $(FFHASH_FILE) already exists."; \
		msg2="Update it (also refresh license)? [Y/N] "; \
		printf "\n%s %s" "$$msg1" "$$msg2"; \
		read ans; echo ""; \
		if [ "$$ans" = "y" ] || [ "$$ans" = "Y" ]; then \
			echo "Re-downloading $(FFHASH_FILE) ..."; \
			echo "Re-downloading $(LICENSE_FILE) ..."; \
			curl -sSL $(FFHASH_REPO)/$(FFHASH_BRANCH)/$(FFHASH_FILE) \
				-o $(FFHASH_FILE); \
			curl -sSL $(FFHASH_REPO)/$(FFHASH_BRANCH)/LICENSE \
				-o $(LICENSE_FILE); \
		else \
			echo "Keeping existing $(FFHASH_FILE)"; \
			echo "Keeping existing $(LICENSE_FILE)"; \
		fi; \
	else \
		echo "Downloading $(FFHASH_FILE) ..."; \
		echo "Downloading $(LICENSE_FILE) ..."; \
		curl -sSL $(FFHASH_REPO)/$(FFHASH_BRANCH)/$(FFHASH_FILE) \
			-o $(FFHASH_FILE); \
		curl -sSL $(FFHASH_REPO)/$(FFHASH_BRANCH)/LICENSE \
			-o $(LICENSE_FILE); \
	fi; \
	echo ""

clean:
	rm -f $(FFHASH_FILE) $(LICENSE_FILE)
