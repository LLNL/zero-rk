
TARGETS = functionTester randomStateGen

all: $(TARGETS)

.PHONY: functionTester
functionTester:
	cd ./functionTester ; make all


.PHONY: randomStateGen
randomStateGen:
	cd ./randomStateGen ; make all

.PHONY:clean
clean:
	cd ./functionTester ; make clean
	cd ./randomStateGen ; make clean

# .PHONY forces make to forces to make to follow a rule unconditionally, for
# some reason the TARGETS weren't being processed
