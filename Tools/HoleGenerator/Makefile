all:
	@if [ ! -d BIN ]; \
	then \
		mkdir BIN; \
		cd BIN; \
		cmake ../SRC; \
		cd ..; \
    fi
	make -C BIN

clean:
	@rm -rf BIN
	@echo "Clean up the build!!"
