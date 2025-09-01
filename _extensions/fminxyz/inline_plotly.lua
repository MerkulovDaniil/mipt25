-- inline_plotly.lua (Combined for Div and Para - Production)

local M = {}

-- Helper function to read file and create HTML block
local function create_plotly_html_output(html_file_path_str)
  local fileContent
  local file, err_msg = io.open(html_file_path_str, "r")
  if file then
    fileContent = file:read("*all")
    file:close()
    return pandoc.RawBlock("html", fileContent)
  else
    -- Return an error message visible in the document if file read fails
    return pandoc.Para({pandoc.Strong({pandoc.Str("[Plotly Filter Error: Could not read file: " .. html_file_path_str .. ". Error: " .. (err_msg or "unknown") .. "]")})})
  end
end

-- Handler for when the plotly directive is parsed as a Div
function M.Div(el)
  if el.classes then
    for _, class_name in ipairs(el.classes) do
      if class_name == "plotly" then
        if el.content and el.content[1] then
          local htmlFilePath = pandoc.utils.stringify(el.content[1]):gsub("^%s*(.-)%s*$", "%1")
          if htmlFilePath ~= "" and string.match(htmlFilePath, "%.html$") then
             return create_plotly_html_output(htmlFilePath)
          end
        end
        -- If plotly class found but content is not a valid path, return a diagnostic
        return pandoc.Para({pandoc.Emph({pandoc.Str("[Plotly Div found, but content was not a valid .html path or was missing.]")})})
      end
    end
  end
  return nil -- Not a plotly Div we handle, or no 'plotly' class
end

-- Handler for when the plotly directive is parsed as a Para (e.g., in lists)
function M.Para(el)
  if el.content and #el.content >= 3 then -- Min parts: start_directive, filename, end_directive
    local first_part_str = pandoc.utils.stringify(el.content[1])

    if string.match(first_part_str, "^%s*:::{%s*%.plotly%s*}%s*$") then
      local filename_str = ""
      -- Try to extract filename, usually at el.content[2].text or el.content[3].text
      if el.content[2] and el.content[2].t == "Str" then
         filename_str = el.content[2].text
      elseif el.content[3] and el.content[3].t == "Str" then
         filename_str = el.content[3].text
      else
         return nil -- Filename not found in expected structure
      end
      
      filename_str = filename_str:gsub("^%s*(.-)%s*$", "%1") -- Trim whitespace

      local last_part_str = pandoc.utils.stringify(el.content[#el.content])

      if string.match(filename_str, "%.html$") and string.match(last_part_str, "^%s*:::%s*$") then
        return create_plotly_html_output(filename_str)
      end
    end
  end
  return nil -- Not a plotly Para we handle
end

return {
  { Div = M.Div },
  { Para = M.Para }
}