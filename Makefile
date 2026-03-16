# Use system typst if available (required for PDF image support, needs >= 0.14)
SYSTEM_TYPST := $(shell command -v typst 2>/dev/null)
ifdef SYSTEM_TYPST
  export QUARTO_TYPST := $(SYSTEM_TYPST)
endif

.PHONY: render render-site render-lecture

render-site:
	quarto render

render-lecture:
	@test -n "$(L)" || (echo "Usage: make render-lecture L=18" && exit 1)
	quarto render lectures/$(L).qmd
