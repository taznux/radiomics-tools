all:
	@if [ ! -d LIB ]; \
	then \
		mkdir LIB; \
		cd LIB; \
		cmake ../SRC; \
		cd ..; \
    fi
	make -C LIB

clean:
	@rm -rf LIB
	@echo "Clean up the build!!"
