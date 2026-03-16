-- Shortcodes for Touying presentation features in typst format
-- Usage in .qmd: {{< pause >}}, {{< meanwhile >}}
-- Must use RawBlock (not RawInline) so Touying treats pauses as block-level separators

function pause()
  if quarto.doc.is_format("typst") then
    return pandoc.RawBlock("typst", "#pause")
  end
end

function meanwhile()
  if quarto.doc.is_format("typst") then
    return pandoc.RawBlock("typst", "#meanwhile")
  end
end
