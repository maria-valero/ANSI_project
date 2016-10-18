
struct seis_header
   {
	int	nt;
	int	id;
	int	rline;	/* receiver line # */
	int	rpt;	/* index on receiver line */
	char	name[8];
	double	x, y, z;
	double	t0;
	double	dt;
   };
#define SEIS_HEADER_SZ	sizeof(struct seis_header)
