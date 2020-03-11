.PHONY: clean All

All:
	@echo "----------Building project:[ new-hpc - Debug ]----------"
	@cd "new-hpc" && "$(MAKE)" -f  "new-hpc.mk"
clean:
	@echo "----------Cleaning project:[ new-hpc - Debug ]----------"
	@cd "new-hpc" && "$(MAKE)" -f  "new-hpc.mk" clean
