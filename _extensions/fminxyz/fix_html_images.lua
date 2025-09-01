-- This Lua filter modifies image paths for HTML output.
-- It replaces PDF images with SVG when available and fixes paths
-- It is designed to be used with Quarto.

function Image (img)
  if FORMAT == "html5" or FORMAT == "html" then
    local original_src = img.src
    local new_src = original_src
    
    -- Handle relative paths - if it doesn't start with http/https, /, or ../, 
    -- assume it's a file in the files directory
    if not new_src:match("^https?://") and not new_src:match("^/") and not new_src:match("^%.%.?/") then
      new_src = "../files/" .. new_src
    end
    
    -- Try to replace PDF with SVG if available
    if new_src:match("%.pdf$") then
      local svg_path = new_src:gsub("%.pdf$", ".svg")
      
      -- Check if the SVG file exists using the same path as in the HTML
      local svg_file = io.open(svg_path, "r")
      if svg_file then
        svg_file:close()
        new_src = svg_path
      else
        print("SVG not found for " .. original_src .. ", keeping PDF")
      end
    end
    
    img.src = new_src
  end
  return img
end