.PHONY: clean All

All:
	@echo "----------Building project:[ hpc-cw - Debug ]----------"
	@cd "hpc-cw" && "$(MAKE)" -f  "hpc-cw.mk"
clean:
	@echo "----------Cleaning project:[ hpc-cw - Debug ]----------"
	@cd "hpc-cw" && "$(MAKE)" -f  "hpc-cw.mk" clean
