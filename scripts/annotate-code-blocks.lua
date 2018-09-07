local vars = {}

function CodeBlock (elem)
    title = nil

    if elem.identifier then
        if vars[elem.identifier] then
            title = pandoc.Str ("«" .. elem.identifier .. "»=+")
        else
            vars[elem.identifier] = true
            title = pandoc.Str ("«" .. elem.identifier .. "»=")
        end
    end

    for k, v in pairs(elem.attr[3]) do
        if k == "file" then
            title = pandoc.Str ("file: «" .. v .. "»=")
        end
    end
    elem.attr[3] = {}

    return { pandoc.Para {pandoc.Emph (title)},  elem }
end
