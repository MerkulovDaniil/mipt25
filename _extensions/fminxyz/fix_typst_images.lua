-- Fix image paths for typst output.
-- Images without explicit path prefix are assumed to be in ../files/

function Image(img)
  if FORMAT == "typst" then
    local src = img.src
    if not src:match("^https?://") and not src:match("^/") and not src:match("^%.%.?/") then
      img.src = "../files/" .. src
    end
  end
  return img
end
