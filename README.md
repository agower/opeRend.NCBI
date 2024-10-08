# EntityClass definitions in opeRend.NCBI

This R package contains EntityClass definitions that correspond to GEO Platform, Sample, and Series objects (for more information, please refer to NCBI's [GEO Overview](https://www.ncbi.nlm.nih.gov/geo/info/overview.html)).

## GEOPlatform
The GEOPlatform entity class captures GEO Platform (GPL) record fields:
| Name          | Array | Type    | Description |
|---------------|-------|---------|-------------|
| geo_accession | FALSE | Text    | 'Platform_geo_accession' field of GEO record (e.g., 'GPL1234') |
| title         | FALSE | Text    | 'Platform_title' field of GEO record |
| manufacturer  | FALSE | Text    | 'Platform_manufacturer' field of GEO record |
| technology    | FALSE | Text    | 'Platform_technology' field of GEO record |
| organism      | TRUE  | Text    | 'Platform_organism' field(s) of GEO record |
| taxid         | TRUE  | Integer | 'Platform_taxid' field of(s) GEO record |

* The `organism` and `taxid` variables are array-type because a Platform record may contain more than one entry for each of these fields (i.e., a platform may include features from multiple organisms).

## GEOSample
The GEOSample entity class captures GEO Sample (GSM) record fields and optionally refers to a CEL file:
| Name                | Array | Type | Description |
|---------------------|-------|------|-------------|
| geo_accession       | FALSE | Text | 'Sample_geo_accession' field of GEO record (e.g., 'GSM1234') |
| title               | FALSE | Text | 'Sample_title' field of GEO record |
| description         | TRUE  | Text | 'Sample_description' field(s) of GEO record |
| type                | FALSE | Text | 'Sample_type' field of GEO record |
| series_id           | TRUE  | Text | 'Sample_series_id' field(s) of GEO record (e.g., 'GSE1234') |
| platform_id         | FALSE | Text | 'Sample_platform_id' field of GEO record (e.g., 'GPL1234') |
| data_processing     | TRUE  | Text | 'Sample_data_processing' field(s) of GEO record |
| channel_count       | FALSE | Integer | 'Sample_channel_count' field of GEO record |
| molecule_ch1        | FALSE | Text | 'Sample_molecule_ch1' field of GEO record |
| source_name_ch1     | FALSE | Text | 'Sample_source_name_ch1' field of GEO record |
| characteristics_ch1 | TRUE  | Text | 'Sample_characteristics_ch1' field(s) of GEO record |
| molecule_ch2        | FALSE | Text | 'Sample_molecule_ch2' field of GEO record |
| source_name_ch2     | FALSE | Text | 'Sample_source_name_ch2' field of GEO record |
| characteristics_ch2 | TRUE  | Text | 'Sample_characteristics_ch2' field(s) of GEO record |
| affymetrixCEL       | FALSE | Entity (AffymetrixCEL) | ID of an AffymetrixCEL Entity that refers to the CEL WorkFile corresponding to this Sample |

* The `description`, `series_id`, and `data_processing` variables are array-type because a Sample record may contain more than one entry for each of these fields.
* The Sample record may contain information about one or two channels, depending on the platform (the number of channels is contained in variable `channel_count`). The `*_ch2` variables therefore may not be present in a given GEOSample record.
* The `characteristics` variables for each channel are array-type because a Sample record may contain more than one characteristics for each channel.
* This record may refer to an AffymetrixCEL entity that describes a CEL WorkFile.

## GEOSeries
The GEOSample entity class captures GEO Series (GSE) record fields, refers to one or more GEOPlatform and GEOSample entities, and may optionally refer to a set of CEL files:
| Name             | Array | Type | Description |
|------------------|-------|------|-------------|
| geo_accession    | FALSE | Text | 'Series_geo_accession' field of GEO record (e.g., 'GSE1234') |
| title            | FALSE | Text | 'Series_title' field of GEO record |
| relation         | TRUE  | Text | 'Series_relation' field of GEO record |
| summary          | TRUE  | Text | 'Series_summary' field(s) of GEO record |
| type             | TRUE  | Text | 'Series_type' field(s) of GEO record (e.g., 'Expression profiling by array') |
| platform_id      | TRUE  | Text | 'Series_platform_id' field(s) of GEO record (e.g., 'GPL1234') |
| pubmed_id        | TRUE  | Integer | 'Series_pubmed_id' field(s) of GEO record |
| geoPlatforms     | TRUE  | Entity (GEOPlatform) | IDs of GEOPlatform Entities corresponding to this Series |
| geoSamples       | TRUE  | Entity (GEOSample) | IDs of GEOSample Entities corresponding to this Series |
| affymetrixCELSet | TRUE  | Entity (AffymetrixCELSet) | IDs of AffymetrixCELSet Entities corresponding to this Series |

* The `relation`, `summary`, `type`, `platform_id`, and `pubmed_id` variables are array-type because a Series record may contain more than one entry for each of these fields.
* Similarly, the entity may refer to more than one GEOPlatform (e.g., methylation and expression arrays).