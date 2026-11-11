library("move2")
library("sf")
library("dplyr")
library("zip")

## The parameter "data" is reserved for the data object passed on from the previous App.
## Settings (ids as in appspec.json): file.name, pts, seg, lins, gpkg

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Drop all attribute columns that only contain NA (the geometry column is kept).
drop_empty_cols <- function(x) {
  x %>% dplyr::select(where(~ !all(is.na(.))))
}

# Convert a move2/sf object into a plain sf table that GDAL can write:
# - removes the move2 class and its track-data/metadata attributes,
# - drops additional (non-active) geometry columns,
# - converts list columns to comma-separated text (GDAL cannot write them).
to_plain_sf <- function(x) {
  gcol <- attr(x, "sf_column")
  geom <- sf::st_geometry(x)
  cols <- setdiff(names(x), gcol)
  out <- list()
  for (n in cols) {
    v <- x[[n]]
    if (inherits(v, "sfc")) next
    if (is.list(v)) {
      v <- vapply(v, function(e) paste(as.character(unlist(e)), collapse = ","), character(1))
      logger.info(paste0("Attribute '", n, "' is a list column and was converted to text for writing."))
    }
    out[[n]] <- v
  }
  df <- structure(out, class = "data.frame", row.names = .set_row_names(length(geom)))
  sf::st_sf(df, geometry = geom)
}

# Shorten attribute names to max. 10 characters for the shapefile format.
# Movebank attributes are renamed with the vocabulary file; other names longer
# than 10 characters keep their first 10 non-vowel characters. Duplicates are
# shortened to 8 characters and indexed. (Same rules as in v3.0.)
shorten_names <- function(nms, vocab) {
  vocab_name  <- tolower(make.names(vocab$name, allow_ = FALSE))
  vocab_short <- make.names(vocab$shortname_mx10, allow_ = FALSE)
  
  nms <- tolower(make.names(nms, allow_ = FALSE))
  novowel_max10 <- substring(gsub("[aeiou ]", "", nms), 1, 10)
  
  for (i in seq_along(nms)) {
    if (nchar(nms[i]) > 10) {
      ix <- match(nms[i], vocab_name)
      nms[i] <- if (!is.na(ix)) vocab_short[ix] else novowel_max10[i]
    }
  }
  
  # if there are duplicates, shorten and index
  X <- 0
  L <- 0
  while (any(duplicated(nms)) & (L + X) < 100) {
    dpl <- which(duplicated(nms))
    for (k in seq_along(dpl)) {
      dd <- which(nms == nms[dpl[k]])
      L <- length(dd)
      nms[dd] <- paste0(substr(nms[dpl[k]], 1, 8), c(1:L + X))
      X <- X + 1
    }
  }
  nms
}

# Rename the attribute columns of a plain sf table using a name map
# (named character vector: original name -> shapefile name).
rename_for_shp <- function(x, name_map, vocab) {
  gcol <- attr(x, "sf_column")
  cols <- setdiff(names(x), gcol)
  new <- unname(name_map[cols])
  if (any(is.na(new))) new[is.na(new)] <- shorten_names(cols[is.na(new)], vocab)
  names(x)[match(cols, names(x))] <- new
  x
}

# Write a plain sf table as a zipped shapefile artefact (folder <file.name>/ inside the zip).
write_zipped_shp <- function(x, shp_dir, file.name, zip_name) {
  unlink(shp_dir, recursive = TRUE)
  sf::st_write(x, shp_dir, layer = file.name, driver = "ESRI Shapefile", quiet = TRUE)
  zipfile <- appArtifactPath(zip_name)
  if (file.exists(zipfile)) unlink(zipfile)
  zip::zip(zipfile = zipfile, files = shp_dir, mode = "cherry-pick")
}

# Write a plain sf table as a GeoPackage artefact.
write_gpkg <- function(x, gpkg_name) {
  f <- appArtifactPath(gpkg_name)
  if (file.exists(f)) unlink(f)
  sf::st_write(x, f, driver = "GPKG", quiet = TRUE)
}

# ------------------------------------------------------------------------------
# Main App function
# ------------------------------------------------------------------------------

rFunction <- function(data, file.name = "moveapps-shapefile", pts = TRUE, seg = FALSE,
                      lins = FALSE, gpkg = FALSE, ...) {
  
  if (is.null(file.name) || is.na(file.name) || trimws(file.name) == "") {
    file.name <- "moveapps-shapefile"
    logger.info("No file name provided, using 'moveapps-shapefile'.")
  }
  file.name <- trimws(file.name)
  
  if (!(pts | seg | lins)) {
    logger.warn(paste("No output type (points, line segments, track lines) was selected.",
                      "No shapefiles", if (gpkg) "or GeoPackages", "are created.",
                      "The input data are returned unchanged."))
    return(data)
  }
  
  ## all track attributes are added to the event table, all-NA attributes are dropped
  pts_m2 <- mt_as_event_attribute(data, names(mt_track_data(data)))
  pts_m2 <- drop_empty_cols(pts_m2)
  pts_sf <- to_plain_sf(pts_m2)
  
  ## shapefile attribute names (max. 10 characters), one map for all outputs
  vocab <- read.csv(getAuxiliaryFilePath("vocab"), header = TRUE)
  attr_names <- setdiff(names(pts_sf), attr(pts_sf, "sf_column"))
  name_map <- setNames(shorten_names(attr_names, vocab), attr_names)
  
  shp_dir <- file.path(tempdir(), file.name)
  logger.info(paste0("Storing shapefiles in tmp-dir ", shp_dir))
  
  ## points
  if (pts) {
    logger.info("You have requested a shapefile with all data points. They are zipped into a file called data_points_shps.zip.")
    write_zipped_shp(rename_for_shp(pts_sf, name_map, vocab), shp_dir, file.name, "data_points_shps.zip")
    if (gpkg) {
      logger.info("You have requested a gpkg file for the points of the track(s).")
      write_gpkg(pts_sf, paste0(file.name, "_points.gpkg"))
    }
  }
  
  ## segments: line between each pair of successive locations, carrying the
  ## attributes of the first location; the last location of each track has no segment
  if (seg) {
    logger.info("You have requested a shapefile with all track segments. They are zipped into a file called data_segments_shps.zip.")
    segs <- mt_segments(pts_m2)
    seg_sf <- pts_sf
    sf::st_geometry(seg_sf) <- segs
    seg_sf <- seg_sf[sf::st_geometry_type(segs) == "LINESTRING", ]
    if (nrow(seg_sf) == 0) {
      logger.warn("No track has more than one location, so no line segments can be created. Segment files are not written.")
    } else {
      write_zipped_shp(rename_for_shp(seg_sf, name_map, vocab), shp_dir, file.name, "data_segments_shps.zip")
      if (gpkg) {
        logger.info("You have requested a gpkg file for the segments of the track(s).")
        write_gpkg(seg_sf, paste0(file.name, "_segments.gpkg"))
      }
    }
  }
  
  ## track lines: one line per track, with all (not all-NA) track attributes
  if (lins) {
    logger.info("You have requested a shapefile with all tracks as lines. They are zipped into a file called data_lines_shps.zip.")
    n_loc <- table(as.character(mt_track_id(data)))
    if (any(n_loc < 2)) {
      logger.warn(paste("The following track(s) have only one location and result in an invalid line:",
                        paste(names(n_loc)[n_loc < 2], collapse = ", ")))
    }
    lins_sf <- drop_empty_cols(to_plain_sf(mt_track_lines(data)))
    write_zipped_shp(rename_for_shp(lins_sf, name_map, vocab), shp_dir, file.name, "data_lines_shps.zip")
    if (gpkg) {
      logger.info("You have requested a gpkg file for the line(s) of the track(s).")
      write_gpkg(lins_sf, paste0(file.name, "_lines.gpkg"))
    }
  }
  
  unlink(shp_dir, recursive = TRUE)
  
  ## the input data are passed on unchanged
  return(data)
}