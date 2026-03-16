#show: fmin-theme.with(
  aspect-ratio: "16-9",
$if(handout)$
  handout: $handout$,
$endif$
$if(bg-image)$
  bg-image: "$bg-image$",
$endif$
$if(title)$
  title: [$title$],
$endif$
$if(by-author)$
  author: [$for(by-author)$$it.name.literal$$sep$, $endfor$],
$endif$
$if(institute)$
  institution: [$institute$],
$endif$
)

#title-slide()
